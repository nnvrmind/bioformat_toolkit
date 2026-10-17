from typing import Iterator

from matplotlib import pyplot as plt

from readers.fastq_record import FastqRecord
from readers.sequence_reader import SequenceReader


class FastqReader(SequenceReader):
    """
    Класс для чтения файлов формата FASTQ.

    Наследуется от SequenceReader и возвращает записи FastqRecord.

    Parameters
    ----------
    filepath: str | Path
        Путь к FASTQ-файлу.
    """

    def __iter__(self) -> Iterator[FastqRecord]:
        """
        Последовательно читает записи FASTQ-файла.

        Yields
        ------
        FastqRecord
            Объект, содержащий идентификатор, последовательность
            и строку качества.

        Raises
        ------
        ValueError
            Если структура FASTQ-записи некорректна.
        FileNotFoundError
            Если входной файл не существует.
        """
        with self.filepath.open("r", encoding="utf-8") as file:
            record_number = 0

            while True:
                header_line = file.readline()

                if header_line == "":
                    return

                record_number += 1
                header = header_line.rstrip("\r\n")

                if not header.startswith("@") or not header[1:]:
                    raise ValueError(
                        f"Запись {record_number}: "
                        "заголовок должен начинаться с @ "
                        "и содержать идентификатор."
                    )

                sequence_lines: list[str] = []

                while True:
                    line = file.readline()

                    if line == "":
                        raise ValueError(
                            f"Запись {record_number}: "
                            "не найден разделитель +."
                        )

                    line = line.rstrip("\r\n")

                    if line.startswith("+"):
                        break

                    sequence_lines.append(line)

                sequence = "".join(sequence_lines)

                if not sequence:
                    raise ValueError(
                        f"Запись {record_number}: "
                        "последовательность не может быть пустой."
                    )

                quality_lines: list[str] = []
                quality_length = 0

                while quality_length < len(sequence):
                    quality_line = file.readline()

                    if quality_line == "":
                        raise ValueError(
                            f"Запись {record_number}: "
                            "строка качества неполная."
                        )

                    quality_line = quality_line.rstrip("\r\n")
                    quality_lines.append(quality_line)
                    quality_length += len(quality_line)

                quality = "".join(quality_lines)

                if len(quality) != len(sequence):
                    raise ValueError(
                        f"Запись {record_number}: "
                        "длина строки качества должна совпадать "
                        "с длиной последовательности."
                    )

                yield FastqRecord(
                    seq_id=header[1:],
                    seq=sequence,
                    quality=quality,
                )

    def plot_per_base_quality(self) -> None:
        """
        Строит график среднего качества секвенирования по позициям.

        Используется кодировка Phred+33.

        Raises
        ------
        ValueError
            Если файл не содержит последовательностей или содержит
            символы качества вне диапазона Phred+33.
        """
        quality_sums: list[float] = []
        quality_counts: list[int] = []

        for record in self:
            for position, quality_char in enumerate(record.quality):
                ascii_code = ord(quality_char)

                if not 33 <= ascii_code <= 126:
                    raise ValueError(
                        "Символ качества не соответствует кодировке Phred+33."
                    )

                score = ascii_code - 33

                if position == len(quality_sums):
                    quality_sums.append(0.0)
                    quality_counts.append(0)

                quality_sums[position] += score
                quality_counts[position] += 1

        if not quality_sums:
            raise ValueError(
                "Невозможно построить график: "
                "FASTQ-файл не содержит последовательностей."
            )

        average_scores = []
        for i in range(len(quality_sums)):
            average = quality_sums[i] / quality_counts[i]
            average_scores.append(average)

        positions = range(1, len(average_scores) + 1)

        plt.plot(positions, average_scores, marker="o")
        plt.xlabel("Позиция в последовательности")
        plt.ylabel("Среднее качество Phred+33")
        plt.title("Среднее качество секвенирования по позициям")
        plt.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    def plot_per_base_content(self) -> None:
        """
        Строит график содержания нуклеотидов по позициям.

        Для каждой позиции вычисляет процент A, C, G и T
        среди всех последовательностей, содержащих эту позицию.

        Raises
        ------
        ValueError
            Если FASTQ-файл не содержит последовательностей.
        """
        bases = "ACGT"

        counts = {}
        for base in bases:
            counts[base] = []

        position_counts: list[int] = []

        for record in self:
            sequence = record.seq.upper()

            for position, nucleotide in enumerate(sequence):
                if position == len(position_counts):
                    position_counts.append(0)

                    for base in bases:
                        counts[base].append(0)

                position_counts[position] += 1

                if nucleotide in counts:
                    counts[nucleotide][position] += 1

        if not position_counts:
            raise ValueError(
                "Невозможно построить график: "
                "FASTQ-файл не содержит последовательностей."
            )

        positions = range(1, len(position_counts) + 1)

        plt.figure()

        for base in bases:
            percentages: list[float] = []

            for position in range(len(position_counts)):
                percentage = (
                    counts[base][position]
                    / position_counts[position]
                    * 100
                )
                percentages.append(percentage)

            plt.plot(positions, percentages, label=base)

        plt.xlabel("Позиция в последовательности")
        plt.ylabel("Содержание нуклеотида (%)")
        plt.title("Содержание нуклеотидов по позициям")
        plt.legend(title="Нуклеотид")
        plt.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()