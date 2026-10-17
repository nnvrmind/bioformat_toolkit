from pathlib import Path

import pandas as pd
import pytest

from readers.sam_reader import SamReader

DATA_DIR = Path(__file__).parent / "data" / "sam"


@pytest.fixture
def simple_sam() -> SamReader:
    """SAM-файл с обычными выравниваниями, unmapped-записью и заголовком."""
    return SamReader(DATA_DIR / "simple.sam")


@pytest.fixture
def edge_sam() -> SamReader:
    """SAM-файл с CIGAR-операциями S, D, I для проверки границ выравнивания."""
    return SamReader(DATA_DIR / "edge_cases.sam")

def test_get_header(simple_sam):
    """
    Возвращает заголовок файла в виде списка строк.

    Expected
    --------
    Список из 6 строк, каждая начинается с символа ``@``.

    Checks
    ------
    Тип результата — ``list``; длина равна 6; все элементы — строки,
    начинающиеся с ``@``.
    """
    header = simple_sam.get_header()
    assert isinstance(header, list)
    assert len(header) == 6
    assert all(line.startswith("@") for line in header)


def test_get_header_group(simple_sam):
    """
    Возвращает строки заголовка, относящиеся к указанной группе.

    Expected
    --------
    ``@SQ`` — 2 строки; ``@RG`` — 1 строка; ``"SQ"`` без ``@`` нормализуется
    и даёт те же 2 строки; несуществующая группа ``@XX`` — пустой список.

    Checks
    ------
    Длина результата для каждой группы; корректная нормализация префикса;
    пустой список для отсутствующей группы.
    """
    assert len(simple_sam.get_header_group("@SQ")) == 2
    assert len(simple_sam.get_header_group("@RG")) == 1
    assert len(simple_sam.get_header_group("SQ")) == 2
    assert simple_sam.get_header_group("@XX") == []


def test_iter_records(simple_sam):
    """
    Итератор по выравниваниям SAM-файла.

    Expected
    --------
    5 записей; у каждой 11 полей; первая запись — ``read1`` с ``flag=0``,
    ``pos=100``, ``cigar="10M"``.

    Checks
    ------
    Количество записей; набор ключей в каждой записи; значения полей
    первой записи.
    """
    records = list(simple_sam)
    assert len(records) == 5
    assert set(records[0].keys()) == {
        "qname", "flag", "rname", "pos", "mapq", "cigar",
        "rnext", "pnext", "tlen", "seq", "qual",
    }
    assert records[0]["qname"] == "read1"
    assert records[0]["flag"] == 0
    assert records[0]["pos"] == 100
    assert records[0]["cigar"] == "10M"

def test_fetch_overlap(simple_sam):
    """
    Возвращает записи, пересекающиеся с интервалом.

    Expected
    --------
    ``read1`` занимает позиции 100-109; запрос ``[105, 115]`` пересекается
    с ним частично → ``read1`` присутствует в результате.

    Checks
    ------
    ``qname`` присутствует среди возвращённых записей.
    """
    qnames = {r["qname"] for r in simple_sam.fetch("chr1", 105, 115)}
    assert "read1" in qnames


def test_fetch_boundary_touch(simple_sam):
    """
    Возвращает записи, пересекающиеся с интервалом (включительные границы).

    Expected
    --------
    ``read1`` заканчивается на 109; запрос ``[109, 120]`` касается границы →
    ``read1`` присутствует в результате.

    Checks
    ------
    ``qname`` присутствует среди возвращённых записей.
    """
    qnames = {r["qname"] for r in simple_sam.fetch("chr1", 109, 120)}
    assert "read1" in qnames


def test_fetch_outside(simple_sam):
    """
    Возвращает записи, пересекающиеся с интервалом.

    Expected
    --------
    ``read1`` занимает 100-109; запрос ``[110, 200]`` не пересекается → пусто.

    Checks
    ------
    Список результатов пуст.
    """
    assert list(simple_sam.fetch("chr1", 110, 200)) == []


def test_fetch_invalid_args(simple_sam):
    """
    Проверяет некорректные аргументы интервала и хромосомы.

    Expected
    --------
    ``start < 1``, ``end < start``, хромосома ``"*"``, пустая строка →
    каждый случай вызывает ``ValueError``.

    Checks
    ------
    ``pytest.raises(ValueError)`` для каждого некорректного вызова.
    """
    with pytest.raises(ValueError):
        list(simple_sam.fetch("chr1", 0, 100))
    with pytest.raises(ValueError):
        list(simple_sam.fetch("chr1", 100, 50))
    with pytest.raises(ValueError):
        list(simple_sam.fetch("*", 1, 100))
    with pytest.raises(ValueError):
        list(simple_sam.fetch("", 1, 100))


def test_fetch_softclip(edge_sam):
    """
    Проверяет учёт soft clip (S) при определении интервала выравнивания.

    Expected
    --------
    ``3S7M`` на позиции 100 занимает 100-106 (soft clip не потребляет
    референс); запрос ``[100, 106]`` попадает, ``[107, 199]`` — нет.

    Checks
    ------
    ``qname`` присутствует/отсутствует среди результатов.
    """
    assert {r["qname"] for r in edge_sam.fetch("chr1", 100, 106)} == {"softclip_start"}
    assert list(edge_sam.fetch("chr1", 107, 199)) == []


def test_fetch_deletion(edge_sam):
    """
    Проверяет учёт deletion (D) при определении интервала выравнивания.

    Expected
    --------
    ``5M2D5M`` на позиции 300 занимает 300-311 (deletion расширяет интервал);
    запрос ``[311, 320]`` попадает, ``[312, 399]`` — нет.

    Checks
    ------
    ``qname`` присутствует/отсутствует среди результатов.
    """
    assert {r["qname"] for r in edge_sam.fetch("chr1", 311, 320)} == {"with_deletion"}
    qnames = {r["qname"] for r in edge_sam.fetch("chr1", 312, 399)}
    assert "with_deletion" not in qnames


def test_fetch_insertion(edge_sam):
    """
    Проверяет учёт insertion (I) при определении интервала выравнивания.

    Expected
    --------
    ``5M2I5M`` на позиции 400 занимает 400-409 (insertion не расширяет
    референсный интервал); запрос ``[400, 409]`` попадает, ``[410, 499]`` — нет.

    Checks
    ------
    ``qname`` присутствует/отсутствует среди результатов.
    """
    assert {r["qname"] for r in edge_sam.fetch("chr1", 400, 409)} == {"with_insertion"}
    assert list(edge_sam.fetch("chr1", 410, 499)) == []


def test_count_alignments(simple_sam):
    """
    Возвращает общее количество выравниваний в файле.

    Expected
    --------
    ``simple.sam`` содержит 5 выравниваний (включая unmapped-запись).

    Checks
    ------
    Возвращаемое значение равно 5.
    """
    assert simple_sam.count_alignments() == 5


def test_count_by_chromosome(simple_sam):
    """
    Возвращает статистику количества выравниваний по хромосомам.

    Expected
    --------
    ``DataFrame`` с колонками ``chromosome``/``count``; ``chr1`` = 3,
    ``chr2`` = 1, ``*`` = 1; строки отсортированы по убыванию ``count``.

    Checks
    ------
    Тип результата; имена колонок; словарь значений; порядок сортировки.
    """
    df = simple_sam.count_by_chromosome()
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["chromosome", "count"]

    counts = dict(zip(df["chromosome"], df["count"]))
    assert counts["chr1"] == 3
    assert counts["chr2"] == 1
    assert counts["*"] == 1

    values = df["count"].tolist()
    assert values == sorted(values, reverse=True)