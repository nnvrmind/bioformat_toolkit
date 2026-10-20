import argparse
from pathlib import Path
from itertools import islice
import pandas as pd

from readers import FastaReader, FastqReader, SamReader, VcfReader


def main() -> None:
    """
    Запускает демонстрационную программу bioformat_toolkit.

    Обрабатывает аргументы командной строки, определяет формат
    входного файла (FASTA, FASTQ, SAM или VCF) и вызывает
    соответствующую демонстрационную функцию.

    Обрабатывает ошибки ввода и работы с файлами.

    Returns
    -------
    None
        Результаты работы программы выводятся в терминал.
    """
    try:
        parser = argparse.ArgumentParser(description="Демонстрационная программа bioformat_toolkit")

        parser.add_argument(
            "format",
            choices=["fasta", "fastq", "sam", "vcf"],
            help="Формат биологического файла"
        )

        parser.add_argument(
            "filepath",
            type=Path,
            help="Путь к файлу"
        )

        args = parser.parse_args()

        if args.format == "fasta":
            demo_fasta(args.filepath)
        elif args.format == "fastq":
            demo_fastq(args.filepath)
        elif args.format == "sam":
            demo_sam(args.filepath)
        elif args.format == "vcf":
            demo_vcf(args.filepath)

    except FileNotFoundError as error:
        print(f"Файл не найден: {error}")
    except PermissionError as error:
        print(f"Нет доступа к файлу: {error}")
    except ValueError as error:
        print(f"Ошибка данных или ввода: {error}")
    except TypeError as error:
        print(f"Ошибка типа данных: {error}")


def demo_fasta(filepath: Path) -> None:
    """
    Демонстрирует чтение и анализ последовательностей из FASTA-файла.

    Parameters
    ----------
    filepath : Path
        Путь к FASTA-файлу.

    Returns
    -------
    None
        Функция выводит идентификаторы, последовательности, их длины,
        общее количество последовательностей и среднюю длину в терминал.
    """
    fasta_reader = FastaReader(filepath)

    for record in fasta_reader:
        print(f"Идентификатор: {record.seq_id}")
        if len(record.seq) > 100:
            print(f"Первые 100 символов последовательности: {record.seq[:100]}...")
        else:
            print(f"Последовательность: {record.seq}")
        print(f"Длина: {len(record)}")
        print()

    print(f"Количество последовательностей: {fasta_reader.count_sequences()}")
    print(f"Средняя длина последовательностей: {fasta_reader.average_length()}")


def demo_fastq(filepath: Path) -> None:
    """
    Демонстрирует чтение и анализ последовательностей из FASTQ-файла.

    Parameters
    ----------
    filepath : Path
        Путь к FASTQ-файлу.

    Returns
    -------
    None
        Функция выводит информацию о последовательностях,
        их количество, среднюю длину и строит графики.
    """
    fastq_reader: FastqReader = FastqReader(filepath)

    print("\nПервые 10 прочтений:")
    for record in islice(fastq_reader, 10):
        print(f"Идентификатор: {record.seq_id}")
        print(f"Последовательность: {record.seq}")
        print(f"Длина: {len(record)}")
        print(f"Качество: {record.quality}")
        print()

    print(f"Количество последовательностей: {fastq_reader.count_sequences()}")
    print(f"Средняя длина последовательностей: {fastq_reader.average_length()}")

    print("Построение графика качества по позициям...")
    fastq_reader.plot_per_base_quality()

    print("Построение графика нуклеотидного состава...")
    fastq_reader.plot_per_base_content()

    print("Построение распределения длин последовательностей...")
    fastq_reader.plot_length_distribution()


def demo_sam(filepath: Path) -> None:
    """
    Демонстрирует чтение и анализ выравниваний из SAM-файла.

    Parameters
    ----------
    filepath : Path
        Путь к SAM-файлу.

    Returns
    -------
    None
        Функция выводит заголовки, количество выравниваний,
        статистику по хромосомам и записи в заданном интервале.
    """
    sam_reader: SamReader = SamReader(filepath)

    print(f"Заголовки: {sam_reader.get_header()}")

    group: str = input("Введите группу заголовка (например, @SQ): ")
    print(f"Заголовки группы {group}: {sam_reader.get_header_group(group)}")

    print(f"Количество выравниваний в SAM-файле: {sam_reader.count_alignments()}")

    print("Количество выравниваний по референсным последовательностям:")
    stats_df: pd.DataFrame = sam_reader.count_by_chromosome()

    if stats_df.empty:
        print("Выравнивания не найдены.")
    else:
        print(stats_df.to_string(index=False))

    reference: str = input("Введите референсную последовательность (например, chr1 или seq1): ")
    start: int = int(input("Введите начало интервала (от 1): "))
    end: int = int(input("Введите конец интервала (включительно): "))

    alignments: list[dict] = list(sam_reader.fetch(reference, start, end))

    print(f"\nНайдено выравниваний: {len(alignments)}")

    if not alignments:
        print("В заданном интервале выравнивания не найдены.")
    else:
        print("\nПервые 10 выравниваний:")

        alignments_df = pd.DataFrame(alignments[:10])
        print(alignments_df.to_string(index=False))


def demo_vcf(filepath: Path) -> None:
    """
    Демонстрирует основные возможности класса VcfReader.

    Выводит заголовки VCF-файла и информацию по выбранной группе,
    общее количество вариантов, записи в заданном геномном
    интервале и статистику по количеству вариантов.

    Parameters
    ----------
    filepath : Path
        Путь к VCF-файлу.

    Returns
    -------
    None
        Результаты работы методов выводятся в терминал.
    """
    vcf_reader = VcfReader(filepath)
    print(f"Заголовки: {vcf_reader.get_header()}")

    group: str = input("Введите группу заголовка: ")
    print(f"Заголовки группы {group}: {vcf_reader.get_header_group(group)}")

    print(f"Количество вариантов в VCF-файле: {vcf_reader.count_variants()}")

    chromosome: str = input("Введите хромосому (например, chr1): ")
    start: int = int(input("Введите начало интервала (от 1): "))
    end: int = int(input("Введите конец интервала (включительно): "))

    records: list[dict[str, str]] = list(vcf_reader.fetch(chromosome, start, end))
    print(f"\nНайдено вариантов: {len(records)}")

    print(f"\nВарианты в интервале [{start}, {end}]:")

    if not records:
        print("В заданном интервале варианты не найдены.")
    else:
        print("\nПервые 10 вариантов:")
        variants_df: pd.DataFrame = pd.DataFrame(records[:10])
        print(variants_df.to_string(index=False))

    region_df: pd.DataFrame = vcf_reader.count_by_region(chromosome, start, end)

    print("\nСтатистика по количеству вариантов:")
    print(region_df.to_string(index=False))


if __name__ == "__main__":
    main()