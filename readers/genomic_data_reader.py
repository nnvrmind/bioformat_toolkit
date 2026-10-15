from readers.reader import Reader
from abc import abstractmethod
from typing import Iterator

class GenomicDataReader(Reader):
    """
    Абстрактный родительский класс для файлов, содержащих данные,
    связанные с геномными координатами.

    Особенности обработки данных реализуются в классах-наследниках.
    """
    
    @abstractmethod
    def get_header(self) -> list[str]:
        """
        Возвращает заголовок файла в виде списка строк.

        Returns
        -------
        list[str]
            Строки заголовка файла.
        """
        pass

    @abstractmethod
    def get_header_group(self, group: str) -> list[str]:
        """
        Возвращает строки заголовка, относящиеся к указанной группе.

        Parameters
        ----------
        group: str
            Название группы заголовка.

        Returns
        -------
        list[str]
            Строки заголовка, принадлежащие группе.
        """
        pass

    @abstractmethod
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
        """
        pass

