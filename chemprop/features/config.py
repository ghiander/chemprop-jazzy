ADDITIONAL_ATOM_PARAMS = {
    "jazzy": {
        # binning of 0.1 from 0.0 to 1.3
        "sd": [x/10 for x in list(range(14))],

        # binning of 0.1 from 0.0 to 1.3
        "sa": [x/10 for x in list(range(14))],

        # from 0 to 4 lone pairs
        "num_lp": list(range(0, 5))
    },

    "kallisto": {
        # binning of 0.1 from -1.0 to 1.0
        "eeq": [x/10 for x in list(range(-10, 10))],

        # binning of 5 from 0 to 60
        # round to the nearest 5 rule
        "alp": list(range(0, 60, 5))
    }
}
