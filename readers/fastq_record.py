from readers.seq import Seq


class FastqRecord(Seq):
    """
    Класс для хранения одной записи FASTQ.

    Наследуется от Seq и дополнительно хранит информацию
    о качестве секвенирования.

    Parameters
    ----------
    seq_id: str
        Идентификатор последовательности.
    seq: str
        Биологическая последовательность.
    quality: str
        Строка символов качества секвенирования.

    Attributes
    ----------
    quality: str
        Строка символов качества секвенирования.
    """

    def __init__(self, seq_id: str, seq: str, quality: str) -> None:
        """
        Инициализирует запись FASTQ.

        Raises
        ------
        TypeError
            Если quality не является строкой.
        ValueError
            Если длина quality не совпадает с длиной seq.
        """
        if not isinstance(quality, str):
            raise TypeError("Строка качества должна быть строкой.")

        super().__init__(seq_id, seq)

        if len(quality) != len(self.seq):
            raise ValueError(
                "Длина строки качества должна совпадать "
                "с длиной последовательности."
            )

        self.quality: str = quality
