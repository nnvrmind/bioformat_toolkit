"""
Модуль тестирования классов Seq и FastaReader.

Этот файл содержит набор автотестов для проверки корректной работы
класса Seq и чтения FASTA-файлов классом FastaReader.

Запуск тестов:
    python -m pytest tests/test_fasta.py -v
"""

import pytest
from pathlib import Path
from typing import List

from readers.seq import Seq
from readers.fasta_reader import FastaReader


def test_seq_length():
    """
    Проверяет корректность работы метода __len__.

    Expected
    --------
    Метод __len__ возвращает целое число, равное количеству символов
    в биологической последовательности.

    Checks
    ------
    Для последовательности 'ATGC' длина должна быть равна 4.
    """
    seq: Seq = Seq("seq1", "ATGC")
    assert len(seq) == 4, "Длина последовательности 'ATGC' должна быть 4"


def test_seq_str():
    """
    Проверяет корректность работы метода __str__.

    Expected
    --------
    Метод __str__ возвращает строку в формате: идентификатор, 
    перенос строки, последовательность.

    Checks
    ------
    Для Seq("seq1", "ACGT") строковое представление должно быть 'seq1\nACGT'.
    """
    seq: Seq = Seq("seq1", "ACGT")
    expected: str = "seq1\nACGT"
    assert str(seq) == expected, f"Ожидалось '{expected}', получено '{str(seq)}'"


def test_seq_attributes():
    """
    Проверяет сохранение атрибутов при создании объекта.

    Expected
    --------
    Атрибуты seq_id и seq содержат значения, переданные в конструктор.

    Checks
    ------
    При создании Seq("gene1", "ATGC") атрибуты равны 'gene1' и 'ATGC'.
    """
    seq: Seq = Seq("gene1", "ATGC")
    assert seq.seq_id == "gene1", "Атрибут seq_id сохранён неверно"
    assert seq.seq == "ATGC", "Атрибут seq сохранён неверно"


def test_seq_empty_id():
    """
    Проверяет валидацию пустого идентификатора.

    Expected
    --------
    При передаче пустой строки в качестве seq_id выбрасывается ValueError.

    Checks
    ------
    Вызов Seq("", "ACGT") должен вызвать исключение ValueError.
    """
    with pytest.raises(ValueError):
        Seq("", "ACGT")


def test_seq_empty_sequence():
    """
    Проверяет валидацию пустой последовательности.

    Expected
    --------
    При передаче пустой строки в качестве seq выбрасывается ValueError.

    Checks
    ------
    Вызов Seq("seq1", "") должен вызвать исключение ValueError.
    """
    with pytest.raises(ValueError):
        Seq("seq1", "")


def test_seq_wrong_type_id():
    """
    Проверяет валидацию типа идентификатора.

    Expected
    --------
    При передаче не-строки в качестве seq_id выбрасывается TypeError.

    Checks
    ------
    Вызов Seq(123, "ACGT") должен вызвать исключение TypeError.
    """
    with pytest.raises(TypeError):
        Seq(123, "ACGT")


def test_seq_wrong_type_seq():
    """
    Проверяет валидацию типа последовательности.

    Expected
    --------
    При передаче не-строки в качестве seq выбрасывается TypeError.

    Checks
    ------
    Вызов Seq("seq1", 123) должен вызвать исключение TypeError.
    """
    with pytest.raises(TypeError):
        Seq("seq1", 123)


def test_fasta_reader_reads_all_records():
    """
    Проверяет, что FastaReader читает все записи из файла.

    Expected
    --------
    Итерация по объекту FastaReader возвращает все записи,
    присутствующие в FASTA-файле.

    Checks
    ------
    Количество полученных записей из example.fasta должно быть равно 3.
    """
    data_path: Path = Path(__file__).parent / "data" / "example.fasta"
    reader: FastaReader = FastaReader(data_path)
    records: List[Seq] = list(reader)
    assert len(records) == 3, "Должно быть прочитано 3 записи"


def test_fasta_reader_content():
    """
    Проверяет корректность извлечения данных из первой записи.

    Expected
    --------
    Идентификатор и последовательность первой записи соответствуют
    данным в файле example.fasta.

    Checks
    ------
    seq_id первой записи == "seq1", seq == "ACGT".
    """
    data_path: Path = Path(__file__).parent / "data" / "example.fasta"
    reader: FastaReader = FastaReader(data_path)
    records: List[Seq] = list(reader)
    assert records[0].seq_id == "seq1", "Неверный seq_id первой записи"
    assert records[0].seq == "ACGT", "Неверная последовательность первой записи"


def test_fasta_reader_invalid_file(tmp_path):
    """
    Проверяет обработку некорректного FASTA-файла.

    Expected
    --------
    Если файл не содержит символа '>' в начале, FastaReader
    должен выбросить ValueError при итерации.

    Checks
    ------
    При чтении файла с содержимым "ACGT\n" выбрасывается ValueError.
    """
    bad_file: Path = tmp_path / "bad.fasta"
    bad_file.write_text("ACGT\n", encoding="utf-8")
    reader: FastaReader = FastaReader(bad_file)
    
    with pytest.raises(ValueError):
        list(reader)


def test_fasta_reader_count_sequences():
    """
    Проверяет работу унаследованного метода count_sequences().

    Expected
    --------
    Метод возвращает точное количество записей в файле.

    Checks
    ------
    Для example.fasta метод count_sequences() возвращает 3.
    """
    data_path: Path = Path(__file__).parent / "data" / "example.fasta"
    reader: FastaReader = FastaReader(data_path)
    assert reader.count_sequences() == 3, "Метод count_sequences вернул неверное число"


def test_fasta_reader_average_length():
    """
    Проверяет работу унаследованного метода average_length().

    Expected
    --------
    Метод вычисляет среднее арифметическое длин всех последовательностей.

    Checks
    ------
    Для длин 4, 6 и 2 средняя длина должна быть равна 4.0.
    """
    data_path: Path = Path(__file__).parent / "data" / "example.fasta"
    reader: FastaReader = FastaReader(data_path)
    assert reader.average_length() == 4.0, "Метод average_length вернул неверное значение"