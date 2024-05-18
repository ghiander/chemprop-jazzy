def create_choices_list(min_value, max_value, step=1):
    """Includes both min and max in the list."""
    if step == 0:
        raise ValueError('Step cannot be 0')
    if step < 0:
        raise ValueError('Step cannot be lower than 0')
    if step < 1:
        min_value = int(min_value * 10)
        max_value = int(max_value * 10)
        step = int(step * 100)
        return [x/step for x in list(range(min_value, max_value + 1))]
    if step >= 1:
        if not float(step).is_integer():
            raise NotImplementedError('Positive steps must be integers')
        return list(range(min_value, max_value + 1, step))
