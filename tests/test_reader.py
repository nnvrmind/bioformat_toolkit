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
    """Проверяет преобразование строкового пути в Path."""
    reader: DummyReader = DummyReader("example.fasta")

    assert isinstance(reader.filepath, Path)
    assert reader.filepath == Path("example.fasta")

def test_filepath_from_path() -> None:
    """Проверяет сохранение пути, переданного как Path."""
    filepath: Path = Path("example.fasta")
    reader: DummyReader = DummyReader(filepath)

    assert reader.filepath == filepath

def test_reader_is_abstract() -> None:
    """Проверяет невозможность создания экземпляра Reader."""
    with pytest.raises(TypeError):
        Reader("example.fasta")