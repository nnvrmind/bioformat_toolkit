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

DATA_DIR = Path(__file__).parent / "data" / "fasta"


def test_seq_length() -> None:
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


def test_seq_str() -> None:
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


def test_seq_attributes() -> None:
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


def test_seq_empty_id() -> None:
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


def test_seq_empty_sequence() -> None:
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


def test_seq_wrong_type_id() -> None:
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


def test_seq_wrong_type_seq() -> None:
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


def test_fasta_reader_reads_all_records() -> None:
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
    
    reader: FastaReader = FastaReader(DATA_DIR / "example.fasta")
    records: List[Seq] = list(reader)
    assert len(records) == 3, "Должно быть прочитано 3 записи"


def test_fasta_reader_content() -> None:
    """
    Проверяет корректность извлечения данных и склейку многострочных записей.

    Expected
    --------
    Все записи читаются правильно, а последовательность, разбитая
    на несколько строк в файле, склеивается в одну.

    Checks
    ------
    Проверяются id и seq всех трех записей. Для seq1 проверяется склейка.
    """
    reader: FastaReader = FastaReader(DATA_DIR / "example.fasta")
    records: List[Seq] = list(reader)
    
    assert records[0].seq_id == "seq1", "Неверный seq_id первой записи"
    assert records[0].seq == (
        "ACGTATGCCTAGCGATCTGATCGATGCCTACTCATGCAGATCGTACGTCGTACGTGCTGCAT"
        "AGACTATCACTGTTTTCATCTACCGCGCGCGCGCTATCGTAGCATGCATCAGCATCGATTGG"
        "TCAGTG"
    ), "Многострочная последовательность seq1 склеена неверно"
    
    assert records[1].seq_id == "seq2", "Неверный seq_id второй записи"
    assert records[1].seq == "ATGCAT", "Неверная последовательность второй записи"
    
    assert records[2].seq_id == "seq3", "Неверный seq_id третьей записи"
    assert records[2].seq == "GG", "Неверная последовательность третьей записи"


def test_fasta_reader_blank_lines() -> None:
    """
    Проверяет, что пустые строки между записями не ломают чтение.

    Expected
    --------
    FastaReader игнорирует пустые строки и корректно читает записи.

    Checks
    ------
    Из файла blank_lines.fasta читаются 2 записи, первая имеет seq == "ACGT".
    """
    reader: FastaReader = FastaReader(DATA_DIR / "blank_lines.fasta")
    records: List[Seq] = list(reader)
    assert len(records) == 2, "Должно быть прочитано 2 записи"
    assert records[0].seq == "ACGT", "Неверная последовательность в первой записи"


def test_fasta_reader_invalid_header() -> None:
    """
    Проверяет обработку файла без символа '>'.

    Expected
    --------
    Если файл не содержит символа '>' в начале, FastaReader
    должен выбросить ValueError при итерации.

    Checks
    ------
    При чтении invalid_header.fasta выбрасывается ValueError.
    """
    reader: FastaReader = FastaReader(DATA_DIR / "invalid_header.fasta")
    with pytest.raises(ValueError):
        list(reader)


def test_fasta_reader_header_only() -> None:
    """
    Проверяет обработку файла, содержащего только заголовок.

    Expected
    --------
    Если у заголовка нет последовательности, FastaReader
    должен выбросить ValueError.

    Checks
    ------
    При чтении header_only.fasta выбрасывается ValueError.
    """
    reader: FastaReader = FastaReader(DATA_DIR / "header_only.fasta")
    with pytest.raises(ValueError):
        list(reader)


def test_fasta_reader_empty_file() -> None:
    """
    Проверяет обработку пустого файла.

    Expected
    --------
    Если файл пуст, FastaReader должен выбросить ValueError.

    Checks
    ------
    При чтении empty.fasta выбрасывается ValueError.
    """
    reader: FastaReader = FastaReader(DATA_DIR / "empty.fasta")
    with pytest.raises(ValueError):
        list(reader)


def test_fasta_reader_count_sequences() -> None:
    """
    Проверяет работу унаследованного метода count_sequences().

    Expected
    --------
    Метод возвращает точное количество записей в файле.

    Checks
    ------
    Для example.fasta метод count_sequences() возвращает 3.
    """
    reader: FastaReader = FastaReader(DATA_DIR / "example.fasta")
    assert reader.count_sequences() == 3, "Метод count_sequences вернул неверное число"


def test_fasta_reader_average_length() -> None:
    """
    Проверяет работу унаследованного метода average_length().

    Expected
    --------
    Метод вычисляет среднее арифметическое длин всех последовательностей.

    Checks
    ------
    Для длин 126, 6 и 2 средняя длина должна быть равна 44.67.
    """
    reader: FastaReader = FastaReader(DATA_DIR / "example.fasta")
    assert round(reader.average_length(), 2) == 46.0, "Метод average_length вернул неверное значение"