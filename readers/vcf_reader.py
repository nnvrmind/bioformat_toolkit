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

        Raises
        ------
        TypeError
            Если название группы не является строкой.
        ValueError
            Если название группы пустое.
        """
        if not isinstance(group, str):
            raise TypeError("Название группы должно быть строкой.")

        if not group:
            raise ValueError("Название группы не может быть пустым.")

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

        Raises
        ------
        ValueError
            Если отсутствует заголовок #CHROM или количество
            полей записи не соответствует количеству столбцов.
        """
        columns: list[str] = []

        with open(self.filepath) as file:
            for line in file:
                if not line.strip():
                    continue

                if line.startswith("#CHROM"):
                    columns = line.rstrip("\r\n").split("\t")
                    columns[0] = columns[0].lstrip("#")

                elif not line.startswith("#"):
                    if not columns:
                        raise ValueError("В VCF-файле отсутствует заголовок #CHROM.")

                    values: list[str] = line.rstrip("\r\n").split("\t")

                    if len(columns) != len(values):
                        raise ValueError("Количество полей записи не соответствует заголовку VCF.")

                    yield dict(zip(columns, values))

        if not columns:  # если файл закончился, а заголовок так и не встретился
            raise ValueError("В VCF-файле отсутствует заголовок #CHROM.")

    def count_variants(self) -> int:
        """
        Возвращает количество вариантов в VCF-файле.

        Returns
        -------
        int
            Общее количество вариантов в файле.
        """
        return sum(1 for _ in self)

    def fetch(self, chromosome: str, start: int, end: int) -> Iterator[dict[str, str]]:
        """
        Последовательно возвращает варианты, позиции которых входят в указанный геномный интервал.

        Parameters
        ----------
        chromosome: str
            Имя хромосомы.
        start: int
            Начало интервала (1-based, включительно).
        end: int
            Конец интервала (1-based, включительно).

        Yields
        ------
        dict[str, str]
            Варианты, позиции которых входят в указанный геномный интервал.

        Raises
        ------
        TypeError
            Если имя хромосомы или координаты имеют некорректный тип.
        ValueError
            Если имя хромосомы пустое или границы интервала некорректны.
        """
        if not isinstance(chromosome, str):
            raise TypeError("Имя хромосомы должно быть строкой.")

        if type(start) is not int or type(end) is not int:
            raise TypeError("Координаты интервала должны быть целыми числами.")

        if not chromosome:
            raise ValueError("Имя хромосомы не может быть пустым.")

        if start < 1 or end < start:
            raise ValueError("Некорректные границы геномного интервала.")

        for record in self:
            if chromosome == record["CHROM"] and start <= int(record["POS"]) <= end:
                yield record

    def count_by_region(self, chromosome: str, start: int, end: int) -> int:
        """
        Возвращает количество вариантов в заданном геномном интервале.

        Parameters
        ----------
        chromosome: str
            Имя хромосомы.
        start: int
            Начало интервала (1-based, включительно).
        end: int
            Конец интервала (1-based, включительно).

        Returns
        -------
        int
            Количество вариантов в указанном интервале.

        Raises
        ------
        TypeError
            Если имя хромосомы или координаты имеют некорректный тип.
        ValueError
            Если имя хромосомы пустое, границы интервала некорректны
            или VCF-файл имеет некорректную структуру.
        """
        records: list[dict[str, str]] = list(self.fetch(chromosome, start, end))
        df: pd.DataFrame = pd.DataFrame(records)
        return len(df)






