from pathlib import Path

from readers.vcf_reader import VcfReader


def test_iter_first_record() -> None:
    """Проверяет чтение первой записи VCF-файла."""
    filepath: Path = Path("tests/data/example.vcf")
    reader: VcfReader = VcfReader(filepath)

    records: list[dict[str, str]] = list(reader) # вызов reader.__iter__() и сбор всех словарей в список

    assert records[0] == {
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
    }
