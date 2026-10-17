"""
Модуль тестирования классов Seq и FastaReader.

Этот файл содержит набор автотестов для проверки корректной работы
класса Seq и чтения FASTA-файлов классом FastaReader.

Запуск тестов:
    python -m pytest tests/test_fasta.py -v
"""

import pytest
from readers.seq import Seq


def test_seq_length():
    """Проверяет, что метод __len__ возвращает правильную длину."""
    seq = Seq("seq1", "ATGC")
    assert len(seq) == 4, "Длина последовательности 'ATGC' должна быть 4"


def test_seq_str():
    """Проверяет, что метод __str__ возвращает строку в формате FASTA."""
    seq = Seq("seq1", "ACGT")
    expected = "seq1\nACGT"
    assert str(seq) == expected, f"Ожидалось '{expected}', получено '{str(seq)}'"


def test_seq_attributes():
    """Проверяет, что атрибуты seq_id и seq сохраняются корректно."""
    seq = Seq("gene1", "ATGC")
    assert seq.seq_id == "gene1", "Атрибут seq_id сохранён неверно"
    assert seq.seq == "ATGC", "Атрибут seq сохранён неверно"


def test_seq_empty_id():
    """Проверяет, что пустой идентификатор вызывает ValueError."""
    with pytest.raises(ValueError):
        Seq("", "ACGT")


def test_seq_empty_sequence():
    """Проверяет, что пустая последовательность вызывает ValueError."""
    with pytest.raises(ValueError):
        Seq("seq1", "")


def test_seq_wrong_type_id():
    """Проверяет, что неверный тип seq_id вызывает TypeError."""
    with pytest.raises(TypeError):
        Seq(123, "ACGT")


def test_seq_wrong_type_seq():
    """Проверяет, что неверный тип seq вызывает TypeError."""
    with pytest.raises(TypeError):
        Seq("seq1", 123)


from readers.fasta_reader import FastaReader
from pathlib import Path

DATA_PATH = Path(__file__).parent / "data" / "example.fasta"


def test_fasta_reader_reads_all_records():
    """Проверяет, что FastaReader прочитал все 3 записи из файла."""
    reader = FastaReader(DATA_PATH)
    records = list(reader)
    assert len(records) == 3


def test_fasta_reader_content():
    """Проверяет, что данные первой записи прочитаны верно."""
    reader = FastaReader(DATA_PATH)
    records = list(reader)
    assert records[0].seq_id == "seq1"
    assert records[0].seq == "ACGT"


def test_fasta_reader_invalid_file(tmp_path):
    """Проверяет, что файл без '>' вызывает ValueError."""
    bad_file = tmp_path / "bad.fasta"
    bad_file.write_text("ACGT\n", encoding="utf-8")
    reader = FastaReader(bad_file)
    
    with pytest.raises(ValueError):
        list(reader)


def test_fasta_reader_count_sequences():
    """Проверяет унаследованный метод count_sequences()."""
    reader = FastaReader(DATA_PATH)
    assert reader.count_sequences() == 3


def test_fasta_reader_average_length():
    """Проверяет унаследованный метод average_length()."""
    reader = FastaReader(DATA_PATH)
    assert reader.average_length() == 4.0