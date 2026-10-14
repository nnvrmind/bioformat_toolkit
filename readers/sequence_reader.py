from readers.reader import Reader


class SequenceReader(Reader):
    """
    Абстрактный промежуточный класс для чтения файлов с последовательностями.

    Наследуется от Reader и добавляет общие для всех форматов
    последовательностей методы: подсчёт количества записей и
    вычисление средней длины.
    """

    def count_sequences(self) -> int:
        """
        Возвращает количество последовательностей в файле.

        Returns
        -------
        int
            Число записей в файле.
        """
        return sum(1 for _ in self)

    def average_length(self) -> float:
        """
        Вычисляет среднюю длину последовательностей в файле.

        Returns
        -------
        float
            Средняя длина последовательностей в файле.
        """
        total_length: int = 0
        count: int = 0

        for record in self:
            total_length += len(record)
            count += 1

        if count == 0:
            return 0.0

        return total_length / count

