from chemprop.features.helpers import create_choices_list


class JazzyValues:
    # binning of 0.1 from 0.0 to 1.3
    sa = {
        "min": 0.0,
        "max": 1.3,
        "step": 0.1
        }
    # binning of 0.1 from 0.0 to 1.3
    sd = {
        "min": 0.0,
        "max": 1.3,
        "step": 0.1
        }
    # from 0 to 4 lone pairs
    num_lp = {
        "min": 0,
        "max": 4,
        "step": 1
        }


class KallistoValues:
    # binning of 0.1 from -1.0 to 0.9
    eeq = {
        "min": -1.0,
        "max": 0.9,
        "step": 0.1
        }
    # binning of 5 from 0 to 55
    # round to the nearest 5 rule
    alp = {
        "min": 0,
        "max": 55,
        "step": 5
        }


ADDITIONAL_ATOM_PARAMS = {
    "jazzy": {
        "sd": create_choices_list(JazzyValues.sd["min"], JazzyValues.sd["max"], JazzyValues.sd["step"]),
        "sa": create_choices_list(JazzyValues.sa["min"], JazzyValues.sa["max"], JazzyValues.sa["step"]),
        "num_lp": create_choices_list(JazzyValues.num_lp["min"], JazzyValues.num_lp["max"], JazzyValues.num_lp["step"])
    },
    "kallisto": {
        "eeq": create_choices_list(KallistoValues.eeq["min"], KallistoValues.eeq["max"], KallistoValues.eeq["step"]),
        "alp": create_choices_list(KallistoValues.alp["min"], KallistoValues.alp["max"], KallistoValues.alp["step"])
    }
}
