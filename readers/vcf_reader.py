from readers.genomic_data_reader import GenomicDataReader
from typing import Iterator
import pandas as pd

class VcfReader(GenomicDataReader):
    """Класс для чтения и анализа файлов формата VCF."""

    def get_header(self) -> list[str]:
        """
        Получает полный заголовок VCF-файла.

        Returns
        -------
        list[str]
            Строки заголовка, включая метаданные и строку с названиями колонок.
        """
        header: list[str] = []
        with open(self.filepath) as file:
            for line in file:
                if line.startswith("#"):
                    header.append(line.rstrip("\r\n"))
                else:
                    break
        return header

    def get_header_group(self, group: str) -> list[str]:
        """
        Получает строки указанной группы метаданных VCF-файла.

        Parameters
        ----------
        group: str
             Название группы метаданных без префикса '##',
             например 'INFO', 'FORMAT' или 'FILTER'.

        Returns
        -------
        list[str]
            Список строк указанной группы.
            Если группа отсутствует, возвращается пустой список.
        """
        header_group: list[str] = []
        with open(self.filepath) as file:
            for line in file:
                if not line.startswith("#"):
                    break
                if line.startswith(f"##{group}="):
                    header_group.append(line.rstrip("\r\n"))
        return header_group

    def __iter__(self) -> Iterator[dict[str, str]]:
        """
        Последовательно читает записи о генетических вариантах.

        Yields
        ------
        dict[str, str]
            Словарь с названиями полей VCF в качестве ключей
            и соответствующими значениями записи.
        """
        with open(self.filepath) as file:
            for line in file:
                if line.startswith("#CHROM"):
                    columns: list[str] = line.rstrip("\r\n").split("\t")
                    columns[0] = columns[0].lstrip("#")
                elif not line.startswith("#"):
                    values: list[str] = line.rstrip("\r\n").split("\t")
                    yield dict(zip(columns, values))

    def count_variants(self) -> int:
        """
        Возвращает количество вариантов в VCF-файле.

        Returns
        -------
        int
            Общее количество вариантов в файле.
        """
        return sum(1 for _ in self)

    def fetch(self, chrom: str, start: int, end: int) -> Iterator[dict[str, str]]:
        """
        Последовательно возвращает варианты, позиции которых входят в указанный геномный интервал.

        Parameters
        ----------
        chrom: str
            Название хромосомы.
        start: int
            Левая граница интервала (включительно).
        end: int
            Правая граница интервала (включительно).

        Yields
        ------
        dict[str, str]
            Варианты, позиции которых входят в указанный геномный интервал.
        """
        for record in self:
            if chrom == record["CHROM"] and start <= int(record["POS"]) <= end:
                yield record

    def count_by_region(self) -> dict[str, int]:
        """
        Подсчитывает количество генетических вариантов на каждой хромосоме.

        Returns
        -------
        dict[str, int]
             Словарь, где ключи — названия хромосом,
             а значения — количество вариантов на них.
        """
        df: pd.DataFrame = pd.DataFrame(self)

        if df.empty:
            return {}

        result: dict[str, int] = df.groupby(by="CHROM").size().to_dict()
        return result






