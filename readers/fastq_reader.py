from typing import Iterator

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