"""
Модуль для чтения биологических последовательностей из файлов формата FASTA.
"""

from pathlib import Path
from typing import Iterator, Union, Optional, List

from readers.seq import Seq
from readers.sequence_reader import SequenceReader


class FastaReader(SequenceReader):
    """
    Класс для чтения FASTA-файлов.
    Наследуется от SequenceReader, получая методы count_sequences и average_length.
    """

    def __iter__(self) -> Iterator[Seq]:
        """
        Читает файл по записям и возвращает объекты Seq.
        Использует генератор (yield), чтобы не загружать весь файл в память.

        Yields
        ------
        Seq
            Очередная запись из FASTA-файла.

        Raises
        ------
        ValueError
            Если файл не соответствует формату FASTA (нет символа '>' в начале).
        """
        with open(self.filepath, "r", encoding="utf-8") as file:
            current_id: Optional[str] = None
            current_seq: List[str] = []

            for line in file:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    if current_id is not None:
                        yield Seq(current_id, "".join(current_seq))
                    
                    current_id = line[1:]
                    current_seq = []
                else:
                    if current_id is None:
                        raise ValueError("Файл не соответствует формату FASTA: нет заголовка '>'")
                    current_seq.append(line)

            if current_id is not None:
                yield Seq(current_id, "".join(current_seq))