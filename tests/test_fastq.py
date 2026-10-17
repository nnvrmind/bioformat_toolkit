import inspect
from pathlib import Path

import pytest
from matplotlib import pyplot as plt

from readers.seq import Seq
from readers.fastq_record import FastqRecord
from readers.fastq_reader import FastqReader


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

FASTQ_PATH = Path(__file__).parent / "data" / "fastq" / "example.fastq"


def test_fastq_reader_reads_three_records() -> None:
    """
    Проверяет чтение всех записей FASTQ-файла.

    Expected
    --------
    Возвращаются три объекта FastqRecord в исходном порядке.

    Checks
    ------
    Проверяются количество записей, идентификаторы и тип объектов.
    """
    reader = FastqReader(FASTQ_PATH)
    records = list(reader)

    assert len(records) == 3
    assert [record.seq_id for record in records] == [
        "EAS54_6_R1_2_1_413_324",
        "EAS54_6_R1_2_1_540_792",
        "EAS54_6_R1_2_1_443_348",
    ]
    assert all(isinstance(record, FastqRecord) for record in records)


def test_fastq_reader_preserves_sequence_and_quality() -> None:
    """
    Проверяет сохранение последовательностей и строк качества.

    Expected
    --------
    Последовательности и качество совпадают с исходным файлом.

    Checks
    ------
    Сравниваются поля полученных объектов с данными FASTQ-файла.
    """
    lines = FASTQ_PATH.read_text(encoding="utf-8").splitlines()

    expected = [
        (lines[i][1:], lines[i + 1], lines[i + 3])
        for i in range(0, len(lines), 4)
    ]

    records = list(FastqReader(FASTQ_PATH))
    actual = [
        (record.seq_id, record.seq, record.quality)
        for record in records
    ]

    assert actual == expected
    assert [len(record) for record in records] == [25, 25, 25]
    assert all(
        len(record.seq) == len(record.quality)
        for record in records
    )


def test_fastq_reader_counts_sequences_and_average_length() -> None:
    """
    Проверяет наследуемые методы подсчёта последовательностей.

    Expected
    --------
    Количество записей равно 3, средняя длина равна 25.0.

    Checks
    ------
    Проверяются результаты count_sequences() и average_length().
    """
    reader = FastqReader(FASTQ_PATH)

    assert reader.count_sequences() == 3
    assert reader.average_length() == 25.0


def test_fastq_reader_returns_generator() -> None:
    """
    Проверяет последовательный характер чтения FASTQ.

    Expected
    --------
    Метод __iter__() возвращает генератор.

    Checks
    ------
    Используется inspect.isgenerator() для проверки результата.
    """
    reader = FastqReader(FASTQ_PATH)

    assert inspect.isgenerator(iter(reader))


def test_fastq_reader_can_be_read_more_than_once() -> None:
    """
    Проверяет повторное чтение файла одним объектом Reader.

    Expected
    --------
    Повторная итерация возвращает те же записи в том же порядке.

    Checks
    ------
    Сравниваются идентификаторы, последовательности и качество.
    """
    reader = FastqReader(FASTQ_PATH)

    first_pass = [
        (record.seq_id, record.seq, record.quality)
        for record in reader
    ]
    second_pass = [
        (record.seq_id, record.seq, record.quality)
        for record in reader
    ]

    assert first_pass == second_pass


@pytest.mark.parametrize(
    "content",
    [
        "read_1\nACGT\n+\nIIII\n",
        "@read_1\nACGT\n-\nIIII\n",
        "@read_1\nACGT\n+\nIII\n",
        "@read_1\nACGT\n+\n",
    ],
)
def test_fastq_reader_rejects_invalid_records(
    tmp_path: Path,
    content: str,
) -> None:
    """
    Проверяет обработку FASTQ-записей с некорректной структурой.

    Expected
    --------
    При неправильной структуре записи возникает ValueError.

    Checks
    ------
    Проверяются неверный заголовок, разделитель, длина качества
    и неполная запись.
    """
    filepath = tmp_path / "invalid.fastq"
    filepath.write_text(content, encoding="utf-8")

    reader = FastqReader(filepath)

    with pytest.raises(ValueError):
        list(reader)


def test_fastq_reader_handles_empty_file(tmp_path: Path) -> None:
    """
    Проверяет чтение пустого FASTQ-файла.

    Expected
    --------
    Итерация не возвращает записей, количество последовательностей равно 0.

    Checks
    ------
    Проверяются list(reader) и count_sequences().
    """
    filepath = tmp_path / "empty.fastq"
    filepath.write_text("", encoding="utf-8")

    reader = FastqReader(filepath)

    assert list(reader) == []
    assert reader.count_sequences() == 0


def test_fastq_reader_raises_for_missing_file(tmp_path: Path) -> None:
    """
    Проверяет обработку отсутствующего файла.

    Expected
    --------
    При попытке чтения несуществующего файла возникает FileNotFoundError.

    Checks
    ------
    Проверяется исключение при запуске итерации.
    """
    filepath = tmp_path / "missing.fastq"
    reader = FastqReader(filepath)

    with pytest.raises(FileNotFoundError):
        list(reader)


def test_fastq_reader_plots_per_base_quality() -> None:
    """
    Проверяет построение графика качества секвенирования по позициям.

    Expected
    --------
    Строится график среднего значения Phred для каждой позиции
    прочтений из FASTQ-файла.

    Checks
    ------
    Проверяется успешное построение графика на тестовом FASTQ-файле.
    """
    plt.show = lambda: None

    reader = FastqReader(FASTQ_PATH)

    reader.plot_per_base_quality()


def test_fastq_reader_plot_empty_file_raises(tmp_path) -> None:
    """
    Проверяет обработку пустого FASTQ-файла при построении графика качества.

    Expected
    --------
    При попытке построить график для пустого файла возникает ValueError.

    Checks
    ------
    Создаётся пустой FASTQ-файл и проверяется возникновение ValueError.
    """
    plt.show = lambda: None

    filepath = tmp_path / "empty.fastq"
    filepath.write_text("", encoding="utf-8")

    reader = FastqReader(filepath)

    with pytest.raises(ValueError):
        reader.plot_per_base_quality()


def test_fastq_reader_plots_per_base_content() -> None:
    """
    Проверяет построение графика содержания нуклеотидов по позициям.

    Expected
    --------
    Строится график с отдельной линией для каждого нуклеотида A, C, G и T.

    Checks
    ------
    Проверяется успешное построение графика на тестовом FASTQ-файле.
    """
    plt.show = lambda: None

    reader = FastqReader(FASTQ_PATH)

    reader.plot_per_base_content()


def test_fastq_reader_content_plot_empty_file_raises(tmp_path) -> None:
    """
    Проверяет обработку пустого FASTQ-файла при построении графика состава.

    Expected
    --------
    При попытке построить график для пустого файла возникает ValueError.

    Checks
    ------
    Создаётся пустой FASTQ-файл и проверяется возникновение ValueError.
    """
    plt.show = lambda: None

    filepath = tmp_path / "empty.fastq"
    filepath.write_text("", encoding="utf-8")

    reader = FastqReader(filepath)

    with pytest.raises(ValueError):
        reader.plot_per_base_content()


def test_fastq_reader_plots_length_distribution() -> None:
    """
    Проверяет построение распределения длин последовательностей.

    Expected
    --------
    Строится график, показывающий количество прочтений каждой длины.

    Checks
    ------
    Проверяется успешное построение графика на тестовом FASTQ-файле.
    """
    plt.show = lambda: None

    reader = FastqReader(FASTQ_PATH)

    reader.plot_length_distribution()


def test_fastq_reader_length_plot_empty_file_raises(tmp_path) -> None:
    """
    Проверяет обработку пустого FASTQ-файла при построении распределения длин.

    Expected
    --------
    При попытке построить график для пустого файла возникает ValueError.

    Checks
    ------
    Создаётся пустой FASTQ-файл и проверяется возникновение ValueError.
    """
    plt.show = lambda: None

    filepath = tmp_path / "empty.fastq"
    filepath.write_text("", encoding="utf-8")

    reader = FastqReader(filepath)

    with pytest.raises(ValueError):
        reader.plot_length_distribution()