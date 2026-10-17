from readers.genomic_data_reader import GenomicDataReader
from typing import Iterator
import pandas as pd
import re

_REF_CONSUMING = frozenset("MDN=X")


def _ref_span(cigar: str) -> int:
    """Длина участка на референсе, занимаемая выравниванием по CIGAR."""
    span = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in _REF_CONSUMING:
            span += int(length)
    return span


class SamReader(GenomicDataReader):
    """
    Чтение SAM-файлов.

    Notes
    -----
    Группы заголовка задаются префиксами ``@HD``, ``@SQ``, ``@RG``,
    ``@PG``, ``@CO``.
    """

    def get_header(self) -> list[str]:
        """
        Возвращает заголовок файла в виде списка строк.

        Returns
        -------
        list[str]
            Строки заголовка файла.
        """
        header: list[str] = []
        with open(self.filepath, "r", encoding="utf-8") as f:
            for line in f:
                if not line.startswith("@"):
                    break
                header.append(line.rstrip("\n"))
        return header

    def __iter__(self) -> Iterator[dict]:
        """
        Итератор по выравниваниям SAM-файла.

        Yields
        ------
        dict
            Словарь с полями выравнивания.
        """
        with open(self.filepath, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("@"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 11:
                    continue
                yield {
                    "qname": fields[0],
                    "flag": int(fields[1]),
                    "rname": fields[2],
                    "pos": int(fields[3]),
                    "mapq": int(fields[4]),
                    "cigar": fields[5],
                    "rnext": fields[6],
                    "pnext": int(fields[7]),
                    "tlen": int(fields[8]),
                    "seq": fields[9],
                    "qual": fields[10],
                }

    def get_header_group(self, group: str) -> list[str]:
        """
        Возвращает строки заголовка, относящиеся к указанной группе.

        Parameters
        ----------
        group: str
            Имя группы заголовка.

        Returns
        -------
        list[str]
            Строки заголовка, принадлежащие группе.
        """
        prefix = group if group.startswith("@") else "@" + group
        return [line for line in self.get_header()
                if line.startswith(prefix)]

    def fetch(self, chromosome: str, start: int, end: int) -> Iterator[dict]:
        """
        Последовательно возвращает записи, расположенные в заданном
        геномном интервале или пересекающиеся с ним.

        Parameters
        ----------
        chromosome: str
            Имя хромосомы (контига).
        start: int
            Начало интервала (1-based, включительно).
        end: int
            Конец интервала (1-based, включительно).

        Yields
        ------
        dict
            Очередная запись, пересекающаяся с интервалом.

        Raises
        ------
        ValueError
            Если интервал задан некорректно.
        """
        if start < 1 or end < start:
            raise ValueError(f"Некорректный интервал: {start}-{end}")
        if not chromosome or chromosome == "*":
            raise ValueError("Хромосома не задана")

        for record in self:
            if record["rname"] != chromosome:
                continue
            if record["cigar"] == "*":
                continue
            if record["flag"] & 0x4:
                continue

            rec_start = record["pos"]                           # 1-based
            rec_end = rec_start + _ref_span(record["cigar"]) - 1  # 1-based, включительно

            if rec_start <= end and rec_end >= start:
                yield record

    def count_alignments(self) -> int:
        """
        Возвращает общее количество выравниваний в файле.

        Returns
        -------
        int
            Число записей-выравниваний (без строк заголовка).
        """
        return sum(1 for _ in self)

    def count_by_chromosome(self) -> pd.DataFrame:
        """
        Возвращает статистику количества выравниваний по хромосомам.

        Returns
        -------
        pandas.DataFrame
            Таблица с колонками ``chromosome`` и ``count``,
            отсортированная по убыванию ``count``.
        """
        chromosomes = [record["rname"] for record in self]
        return (
            pd.Series(chromosomes, name="count")
            .value_counts()
            .rename_axis("chromosome")
            .reset_index(name="count")
        )