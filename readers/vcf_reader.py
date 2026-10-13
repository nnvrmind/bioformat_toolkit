from genomic_data_reader import GenomicDataReader

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






