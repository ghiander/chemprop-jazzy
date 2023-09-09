from chemprop.features.config import ADDITIONAL_ATOM_PARAMS
from chemprop.features.encoding import onek_encoding_unk


def _calculate_param_lengths():
    """Calculates the dimensionalities of the features."""
    length_dict = dict()
    for library, features in ADDITIONAL_ATOM_PARAMS.items():
        length_dict[library] = sum([__len_including_zero(f)
                                    for f in features.values()])
    return length_dict


def __len_including_zero(lst: list):
    return len(lst) + 1


ADDITIONAL_ATOM_PARAM_LENGTHS = _calculate_param_lengths()


def encode_feature_list(library, feat_list):
    """
    One-hot encodes all atomic features of a list
    for a given calculation library (Jazzy).
    """
    for idx in range(len(feat_list)):
        atom_feat_dict = feat_list[idx]
        for feature in atom_feat_dict.keys():
            atom_feat_dict[feature] = _encode_feature(library,
                                                      feature,
                                                      atom_feat_dict)
    return feat_list


def _encode_feature(library, feature, atom_feat_dict):
    return onek_encoding_unk(atom_feat_dict[feature],
                             ADDITIONAL_ATOM_PARAMS[library][feature])
