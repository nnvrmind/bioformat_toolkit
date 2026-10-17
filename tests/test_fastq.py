import pytest

from readers.seq import Seq
from readers.fastq_record import FastqRecord


def test_fastq_record_inherits_seq() -> None:
    """
    Проверяет наследование FastqRecord от Seq.

    Expected
    --------
    FastqRecord является подклассом Seq.

    Checks
    ------
    Проверяется отношение наследования с помощью issubclass().
    """
    assert issubclass(FastqRecord, Seq)


def test_fastq_record_stores_fields() -> None:
    """
    Проверяет создание FastqRecord и сохранение его атрибутов.

    Expected
    --------
    Объект сохраняет идентификатор, последовательность и качество.

    Checks
    ------
    Проверяются значения seq_id, seq и quality.
    """
    record = FastqRecord("read_1", "ACGT", "IIII")

    assert record.seq_id == "read_1"
    assert record.seq == "ACGT"
    assert record.quality == "IIII"


def test_fastq_record_inherits_seq_methods() -> None:
    """
    Проверяет работу унаследованных методов Seq.

    Expected
    --------
    len() возвращает длину последовательности, а str()
    возвращает идентификатор и последовательность на разных строках.

    Checks
    ------
    Проверяются результаты len(record) и str(record).
    """
    record = FastqRecord("read_1", "ACGT", "IIII")

    assert len(record) == 4
    assert str(record) == "read_1\nACGT"


def test_fastq_record_requires_string_quality() -> None:
    """
    Проверяет тип данных атрибута quality.

    Expected
    --------
    При передаче quality нестрокового типа возникает TypeError.

    Checks
    ------
    Передаётся числовое значение вместо строки качества.
    """
    with pytest.raises(TypeError):
        FastqRecord("read_1", "ACGT", 123)


@pytest.mark.parametrize("quality", ["", "III", "IIIII"])
def test_fastq_record_requires_matching_lengths(quality: str) -> None:
    """
    Проверяет соответствие длины quality длине последовательности.

    Expected
    --------
    При несовпадении длин возникает ValueError.

    Checks
    ------
    Проверяются пустая, более короткая и более длинная строки качества.
    """
    with pytest.raises(ValueError):
        FastqRecord("read_1", "ACGT", quality)
