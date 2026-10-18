from pathlib import Path
from typing import Iterator
import pytest
from readers.reader import Reader


class DummyReader(Reader):
    """Вспомогательный класс для тестирования Reader."""

    def __iter__(self) -> Iterator:
        """Возвращает пустой итератор."""
        records: list = []
        return iter(records)


def test_filepath_from_string() -> None:
    """
    Проверяет преобразование строкового пути в объект Path.

    Expected
    --------
    Атрибут filepath содержит объект Path, соответствующий
    переданному строковому пути.

    Checks
    ------
    Тип атрибута filepath и соответствие исходному пути.
    """
    reader: DummyReader = DummyReader("example.fasta")

    assert isinstance(reader.filepath, Path)
    assert reader.filepath == Path("example.fasta")


def test_filepath_from_path() -> None:
    """
    Проверяет сохранение пути, переданного как объект Path.

    Expected
    --------
    Атрибут filepath содержит переданный объект Path.

    Checks
    ------
    Равенство атрибута filepath исходному пути.
    """
    filepath: Path = Path("example.fasta")
    reader: DummyReader = DummyReader(filepath)

    assert reader.filepath == filepath


def test_reader_is_abstract() -> None:
    """
    Проверяет невозможность создания экземпляра абстрактного класса Reader.

    Expected
    --------
    При попытке создать экземпляр Reader возникает TypeError.

    Checks
    ------
    Возникновение исключения TypeError при создании экземпляра.
    """
    with pytest.raises(TypeError):
        Reader("example.fasta")