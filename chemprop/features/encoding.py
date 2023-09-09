from typing import List


def onek_encoding_unk(value: int, choices: List[int]) -> List[int]:
    """
    Creates a one-hot encoding with an extra category for uncommon values.

    :param value: The value for which the encoding should be one.
    :param choices: A list of possible values.
    :return: A one-hot encoding of the :code:`value` in a list of length :code:`len(choices) + 1`.
             If :code:`value` is not in :code:`choices`, then the final element in the encoding is 1.
    """
    encoding = [0] * (len(choices) + 1)
    index = choices.index(value) if value in choices else -1
    encoding[index] = 1

    return encoding


def bin_value(value, binning):
    """
    Bins a value according to a binning factor.
    e.g. binning == 5; 43 becomes 45, 41 becomes 40.
    """
    return binning*round(value/binning)


def bin_fraction(value, binning_ratio):
    """
    Bins a decimal value according to a binning ratio.
    e.g. binning_ratio == 0.05; 0.12 becomes 0.1, 0.03 becomes 0.05.
    """
    binning = binning_ratio*100
    return bin_value(round(value, 2)*100,
                     binning)/100
