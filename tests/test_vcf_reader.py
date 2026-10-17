from pathlib import Path
import pytest
from readers.vcf_reader import VcfReader


def test_get_header() -> None:
    """
    Проверяет корректность получения заголовка VCF-файла.

    Expected
    --------
    Метод get_header() возвращает список из пяти строк заголовка
    в исходном порядке, без символов переноса строки.

    Checks
    ------
    Полное соответствие полученного списка ожидаемому:
    содержимое, количество и порядок строк.
    """
    filepath: Path = Path("tests/data/vcf/example.vcf")
    reader: VcfReader = VcfReader(filepath)

    header: list[str] = reader.get_header()

    expected: list[str] = [
        '##fileformat=VCFv4.2',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        '##FILTER=<ID=q10,Description="Quality below 10">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1'
    ]

    assert header == expected


def test_get_header_group() -> None:
    """
    Проверяет получение групп метаданных из заголовка VCF-файла.

    Expected
    --------
    Метод get_header_group() возвращает строки указанных групп
    метаданных без символов переноса строки.

    Checks
    ------
    Полное соответствие объединённого списка ожидаемому:
    содержимое, количество и порядок строк.
    """
    filepath: Path = Path("tests/data/vcf/example.vcf")
    reader: VcfReader = VcfReader(filepath)

    header_group: list[str] = []
    groups: list[str] = ["fileformat", "INFO", "FILTER", "FORMAT"]
    for group in groups:
        header_group.extend(reader.get_header_group(group))

    expected: list[str] = [
        '##fileformat=VCFv4.2',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        '##FILTER=<ID=q10,Description="Quality below 10">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    ]

    assert header_group == expected


def test_get_header_group_invalid_type() -> None:
    """
    Проверяет обработку некорректного типа аргумента в get_header_group().

    Expected
    --------
    Метод get_header_group() вызывает TypeError,
    если название группы метаданных не является строкой.

    Checks
    ------
    Возникновение исключения TypeError при передаче
    целого числа вместо названия группы.
    """
    filepath: Path = Path("tests/data/vcf/example.vcf")
    reader: VcfReader = VcfReader(filepath)

    with pytest.raises(TypeError):
        reader.get_header_group(123)


def test_get_header_group_missing_group() -> None:
    """
    Проверяет получение отсутствующей группы метаданных VCF-файла.

    Expected
    --------
    Метод get_header_group() возвращает пустой список,
    если указанная группа отсутствует в заголовке файла.

    Checks
    ------
    Соответствие результата пустому списку при запросе
    несуществующей группы FILE.
    """
    filepath: Path = Path("tests/data/vcf/example.vcf")
    reader: VcfReader = VcfReader(filepath)

    header_group: list[str] = reader.get_header_group("FILE")

    assert header_group == []


def test_get_header_group_empty_name() -> None:
    """
    Проверяет обработку пустого названия группы метаданных.

    Expected
    --------
    Метод get_header_group() вызывает ValueError
    при передаче пустой строки.

    Checks
    ------
    Возникновение исключения ValueError.
    """
    filepath: Path = Path("tests/data/vcf/example.vcf")
    reader: VcfReader = VcfReader(filepath)

    with pytest.raises(ValueError):
        reader.get_header_group("")


def test_iter_records() -> None:
    """
    Проверяет корректность чтения всех записей VCF-файла.

    Expected
    --------
    Метод __iter__() возвращает восемь записей VCF-файла
    в виде словарей, где ключи соответствуют названиям
    столбцов, а значения представлены строками.

    Checks
    ------
    Полное соответствие полученного списка ожидаемому:
    количество, содержимое и порядок записей.
    """
    filepath: Path = Path("tests/data/vcf/example.vcf")
    reader: VcfReader = VcfReader(filepath)

    records: list[dict[str, str]] = list(reader)  # вызов reader.__iter__() и сбор всех словарей в список

    expected: list[dict[str, str]] = [
        {
            "CHROM": "chr1",
            "POS": "1000",
            "ID": "rs123",
            "REF": "A",
            "ALT": "G",
            "QUAL": "99",
            "FILTER": "PASS",
            "INFO": "DP=35",
            "FORMAT": "GT",
            "SAMPLE1": "0/1",
        },
        {
            "CHROM": "chr1",
            "POS": "1200",
            "ID": ".",
            "REF": "G",
            "ALT": "C",
            "QUAL": "75",
            "FILTER": "PASS",
            "INFO": "DP=28",
            "FORMAT": "GT",
            "SAMPLE1": "0/1",
        },
        {
            "CHROM": "chr1",
            "POS": "1500",
            "ID": ".",
            "REF": "C",
            "ALT": "T",
            "QUAL": "80",
            "FILTER": "PASS",
            "INFO": "DP=42",
            "FORMAT": "GT",
            "SAMPLE1": "1/1",
        },
        {
            "CHROM": "chr1",
            "POS": "2000",
            "ID": "rs789",
            "REF": "T",
            "ALT": "A",
            "QUAL": "60",
            "FILTER": "PASS",
            "INFO": "DP=31",
            "FORMAT": "GT",
            "SAMPLE1": "0/1",
        },
        {
            "CHROM": "chr2",
            "POS": "1000",
            "ID": ".",
            "REF": "A",
            "ALT": "C",
            "QUAL": "90",
            "FILTER": "PASS",
            "INFO": "DP=50",
            "FORMAT": "GT",
            "SAMPLE1": "0/1",
        },
        {
            "CHROM": "chr2",
            "POS": "2000",
            "ID": "rs456",
            "REF": "G",
            "ALT": "A",
            "QUAL": "50",
            "FILTER": "PASS",
            "INFO": "DP=18",
            "FORMAT": "GT",
            "SAMPLE1": "0/0",
        },
        {
            "CHROM": "chr2",
            "POS": "2500",
            "ID": ".",
            "REF": "C",
            "ALT": "G",
            "QUAL": "8",
            "FILTER": "q10",
            "INFO": "DP=12",
            "FORMAT": "GT",
            "SAMPLE1": "0/1",
        },
        {
            "CHROM": "chr3",
            "POS": "3000",
            "ID": "rs999",
            "REF": "T",
            "ALT": "C",
            "QUAL": "95",
            "FILTER": "PASS",
            "INFO": "DP=40",
            "FORMAT": "GT",
            "SAMPLE1": "1/1",
        },
    ]

    assert records == expected


def test_iter_missing_header() -> None:
    """
    Проверяет обработку отсутствующего заголовка #CHROM в VCF-файле.

    Expected
    --------
    Метод __iter__() вызывает ValueError при отсутствии заголовка #CHROM.

    Checks
    ------
    Возникновение исключения ValueError при чтении VCF-файла
    без заголовка #CHROM.
    """
    filepath: Path = Path("tests/data/vcf/missing_header.vcf")
    reader: VcfReader = VcfReader(filepath)

    with pytest.raises(ValueError, match="отсутствует заголовок"):
        list(reader)


def test_iter_invalid_columns() -> None:
    """
    Проверяет обработку некорректного количества полей записи VCF-файла.

    Expected
    --------
    Метод __iter__() вызывает ValueError, если количество полей записи
    не соответствует количеству столбцов заголовка.

    Checks
    ------
    Возникновение исключения ValueError при чтении VCF-файла
    с некорректным количеством полей.
    """
    filepath: Path = Path("tests/data/vcf/invalid_columns.vcf")
    reader: VcfReader = VcfReader(filepath)

    with pytest.raises(ValueError, match="Количество полей записи"):
        list(reader)


def test_iter_header_only() -> None:
    """
    Проверяет чтение VCF-файла, содержащего только заголовок.

    Expected
    --------
    Метод __iter__() возвращает пустой итератор.

    Checks
    ------
    Отсутствие записей при чтении файла без вариантов.
    """
    filepath: Path = Path("tests/data/vcf/header_only.vcf")
    reader: VcfReader = VcfReader(filepath)

    assert list(reader) == []


def test_iter_blank_lines() -> None:
    """
    Проверяет пропуск пустых строк при чтении VCF-файла.

    Expected
    --------
    Метод __iter__() пропускает пустые строки
    и возвращает только записи генетических вариантов.

    Checks
    ------
    Количество прочитанных записей и их позиции.
    """
    filepath: Path = Path("tests/data/vcf/blank_lines.vcf")
    reader: VcfReader = VcfReader(filepath)

    records: list[dict[str, str]] = list(reader)

    assert len(records) == 2
    assert [record["POS"] for record in records] == ["1000", "2000"]


def test_iter_empty_file() -> None:
    """
    Проверяет обработку пустого VCF-файла.

    Expected
    --------
    При чтении пустого файла возникает ValueError.

    Checks
    ------
    Возникновение исключения ValueError.
    """
    filepath: Path = Path("tests/data/vcf/empty.vcf")
    reader: VcfReader = VcfReader(filepath)

    with pytest.raises(ValueError):
        list(reader)


def test_count_variants() -> None:
    """
    Проверяет подсчёт количества вариантов в VCF-файле.

    Expected
    --------
    Метод count_variants() возвращает 8.

    Checks
    ------
    Соответствие количества вариантов в файле ожидаемому значению.
    """
    filepath: Path = Path("tests/data/vcf/example.vcf")
    reader: VcfReader = VcfReader(filepath)

    assert reader.count_variants() == 8


def test_count_variants_header_only() -> None:
    """
    Проверяет подсчёт вариантов в VCF-файле без записей.

    Expected
    --------
    Метод count_variants() возвращает 0.

    Checks
    ------
    Равенство количества вариантов нулю.
    """
    filepath: Path = Path("tests/data/vcf/header_only.vcf")
    reader: VcfReader = VcfReader(filepath)

    assert reader.count_variants() == 0


def test_fetch() -> None:
    """
    Проверяет получение вариантов из заданного геномного интервала методом fetch().

    Expected
    --------
    При преобразовании результата fetch("chr1", 500, 1000) в список
    получается одна запись VCF с позицией 1000.
    Правая граница интервала включается в результат.

    Checks
    ------
    Полное соответствие полученного списка ожидаемому:
    количество записей, названия полей и их значения.
    """
    filepath: Path = Path("tests/data/vcf/example.vcf")
    reader: VcfReader = VcfReader(filepath)

    records: list[dict[str, str]] = list(reader.fetch("chr1", 500, 1000))
    expected: list[dict[str, str]] = [{
        "CHROM": "chr1",
        "POS": "1000",
        "ID": "rs123",
        "REF": "A",
        "ALT": "G",
        "QUAL": "99",
        "FILTER": "PASS",
        "INFO": "DP=35",
        "FORMAT": "GT",
        "SAMPLE1": "0/1",
    }]
    assert records == expected


def test_fetch_left_boundary() -> None:
    """
    Проверяет включение левой границы геномного интервала в fetch().

    Expected
    --------
    Метод fetch() возвращает один вариант с координатой POS = 1000,
    совпадающей с началом заданного интервала.

    Checks
    ------
    Количество найденных вариантов равно 1.
    Координата найденного варианта равна 1000.
    """
    filepath: Path = Path("tests/data/vcf/example.vcf")
    reader: VcfReader = VcfReader(filepath)

    records: list[dict[str, str]] = list(reader.fetch("chr1", 1000, 1100))

    assert len(records) == 1
    assert records[0]["POS"] == "1000"


def test_fetch_several_records() -> None:
    """
    Проверяет получение нескольких вариантов из заданного геномного интервала.

    Expected
    --------
    При преобразовании результата fetch("chr1", 1000, 2000) в список получается четыре записи VCF,
    расположенные на хромосоме chr1 в позициях 1000, 1200, 1500 и 2000.
    Обе границы интервала включаются в результат.

    Checks
    ------
    Полное соответствие полученного списка ожидаемому:
    количество записей, порядок, названия полей и их значения.
    """
    filepath: Path = Path("tests/data/vcf/example.vcf")
    reader: VcfReader = VcfReader(filepath)

    records: list[dict[str, str]] = list(reader.fetch("chr1", 1000, 2000))
    expected: list[dict[str, str]] = [
        {
            "CHROM": "chr1",
            "POS": "1000",
            "ID": "rs123",
            "REF": "A",
            "ALT": "G",
            "QUAL": "99",
            "FILTER": "PASS",
            "INFO": "DP=35",
            "FORMAT": "GT",
            "SAMPLE1": "0/1",
        },
        {
            "CHROM": "chr1",
            "POS": "1200",
            "ID": ".",
            "REF": "G",
            "ALT": "C",
            "QUAL": "75",
            "FILTER": "PASS",
            "INFO": "DP=28",
            "FORMAT": "GT",
            "SAMPLE1": "0/1",
        },
        {
            "CHROM": "chr1",
            "POS": "1500",
            "ID": ".",
            "REF": "C",
            "ALT": "T",
            "QUAL": "80",
            "FILTER": "PASS",
            "INFO": "DP=42",
            "FORMAT": "GT",
            "SAMPLE1": "1/1",
        },
        {
            "CHROM": "chr1",
            "POS": "2000",
            "ID": "rs789",
            "REF": "T",
            "ALT": "A",
            "QUAL": "60",
            "FILTER": "PASS",
            "INFO": "DP=31",
            "FORMAT": "GT",
            "SAMPLE1": "0/1",
        },
    ]
    assert records == expected


def test_fetch_empty_result() -> None:
    """
    Проверяет получение пустого результата при отсутствии вариантов в интервале.

    Expected
    --------
    При преобразовании результата fetch("chr1", 3000, 4000) в список получается пустой список,
    поскольку в указанном геномном интервале нет вариантов.

    Checks
    ------
    Полученный список не содержит записей VCF.
    """
    filepath: Path = Path("tests/data/vcf/example.vcf")
    reader: VcfReader = VcfReader(filepath)

    records: list[dict[str, str]] = list(reader.fetch("chr1", 3000, 4000))

    assert records == []


# декоратор, позволяющий запускать тестовую функцию несколько раз,
# каждый раз подставляя разные значения её параметров
@pytest.mark.parametrize(
    "chromosome, start, end",
    [
        (123, 1000, 2000),
        ("chr1", "1000", 2000),
        ("chr1", 1000, 2000.5),
        ("chr1", True, 2000),
    ]
)
def test_fetch_invalid_type(
    chromosome: str | int,
    start: int | str | bool,
    end: int | float,
) -> None:
    """
    Проверяет обработку некорректных типов аргументов метода fetch().

    Expected
    --------
    Метод fetch() вызывает TypeError, если имя хромосомы не является
    строкой или координаты интервала не являются целыми числами.

    Checks
    ------
    Возникновение исключения TypeError при передаче:
    - целого числа вместо имени хромосомы;
    - строки вместо начальной координаты;
    - дробного числа вместо конечной координаты;
    - логического значения вместо начальной координаты.
    """
    filepath: Path = Path("tests/data/vcf/example.vcf")
    reader: VcfReader = VcfReader(filepath)

    with pytest.raises(TypeError):
        list(reader.fetch(chromosome, start, end))


@pytest.mark.parametrize(
    "chromosome, start, end",
    [
        ("", 1000, 2000),
        ("chr1", 0, 2000),
        ("chr1", 2000, 1000),
    ]
)
def test_fetch_invalid_value(chromosome: str, start: int, end: int) -> None:
    """
    Проверяет обработку некорректных значений аргументов метода fetch().

    Expected
    --------
    Метод fetch() вызывает ValueError, если имя хромосомы пустое
    или границы геномного интервала некорректны.

    Checks
    ------
    Возникновение исключения ValueError при передаче:
    - пустого имени хромосомы;
    - начальной координаты, равной нулю;
    - начальной координаты, превышающей конечную.
    """
    filepath: Path = Path("tests/data/vcf/example.vcf")
    reader: VcfReader = VcfReader(filepath)

    with pytest.raises(ValueError):
        list(reader.fetch(chromosome, start, end))


def test_count_by_region() -> None:
    """
    Проверяет подсчёт вариантов в заданном геномном интервале.

    Expected
    --------
    Метод count_by_region() возвращает количество вариантов
    на хромосоме chr1 в интервале [1000, 2000], равное 4.

    Checks
    ------
    Соответствие возвращаемого количества ожидаемому значению.
    """
    filepath: Path = Path("tests/data/vcf/example.vcf")
    reader: VcfReader = VcfReader(filepath)

    count: int = reader.count_by_region("chr1", 1000, 2000)
    assert count == 4


def test_count_by_region_no_variants() -> None:
    """
    Проверяет подсчёт вариантов в интервале без записей.

    Expected
    --------
    Метод count_by_region() возвращает 0, если в заданном
    геномном интервале отсутствуют варианты.

    Checks
    ------
    Равенство возвращаемого количества нулю.
    """
    filepath: Path = Path("tests/data/vcf/example.vcf")
    reader: VcfReader = VcfReader(filepath)

    count: int = reader.count_by_region("chr1", 3000, 4000)
    assert count == 0


@pytest.mark.parametrize(
    "chromosome, start, end, expected_exception",
    [
        (123, 1000, 2000, TypeError),
        ("chr1", "1000", 2000, TypeError),
        ("chr1", True, 2000, TypeError),
        ("", 1000, 2000, ValueError),
        ("chr1", 0, 2000, ValueError),
        ("chr1", 2000, 1000, ValueError),
    ]
)
def test_count_by_region_invalid_arguments(
    chromosome: str | int,
    start: int | str | bool,
    end: int,
    expected_exception: type[Exception],
) -> None:
    """
    Проверяет обработку некорректных аргументов в count_by_region().

    Expected
    --------
    Метод count_by_region() вызывает TypeError или ValueError
    при передаче некорректных аргументов.

    Checks
    ------
    Передача исключений из fetch() в count_by_region()
    при неправильных типах и значениях аргументов.
    """
    filepath: Path = Path("tests/data/vcf/example.vcf")
    reader: VcfReader = VcfReader(filepath)

    with pytest.raises(expected_exception):
        reader.count_by_region(chromosome, start, end)