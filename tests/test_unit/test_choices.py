from chemprop.features.helpers import create_choices_list


def test_decimal_choices_list_1():
    choices = create_choices_list(-0.5, 0.5, 0.1)
    expected = [-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
    assert choices == expected


def test_decimal_choices_list_2():
    choices = create_choices_list(0, 1.3, 0.1)
    expected = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
    assert choices == expected


def test_nearest_five_choices_list():
    choices = create_choices_list(0, 60, 5)
    expected = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
    assert choices == expected


def test_choices_list():
    choices = create_choices_list(0, 5)
    expected = [0, 1, 2, 3, 4, 5]
    assert choices == expected
