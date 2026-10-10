from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterator

class Reader(ABC):
    """
    Абстрактный базовый класс для чтения биоинформатических файлов.

    Parameters
    ----------
    filepath: str | Path
        Путь к файлу.

    Attributes
    -----------
    filepath: Path
        Путь к файлу в виде объекта pathlib.Path.
    """

    def __init__(self, filepath: str | Path) -> None:
        """
        Инициализирует объект Reader и сохраняет путь к файлу.

        Parameters
        ----------
        filepath: str | Path
            Путь к файлу.
        """
        self.filepath = Path(filepath)

    @abstractmethod
    def __iter__(self) -> Iterator:
        """
        Возвращает итератор для последовательного чтения записей файла.

        Returns
        -------
        Iterator
            Итератор, выдающий записи файла по одной.
        """
        pass