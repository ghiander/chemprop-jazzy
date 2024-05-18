import pytest
from rdkit import Chem
from chemprop.features.jazzy import calculate_jazzy_and_kallisto_features
from chemprop.features.additional_features import encode_feature_list


def _create_water_mol():
    smi = "[H]O[H]"
    mol = Chem.MolFromSmiles(smi)
    return Chem.AddHs(mol)


def _get_atom_dict(mol):
    return {atom.GetIdx(): atom.GetAtomicNum()
            for atom in mol.GetAtoms()}


def _get_water_oxygen_idx(water_mol):
    atom_dict = _get_atom_dict(water_mol)
    for k, v in atom_dict.items():
        if v == 8:
            return k


def _get_water_hydrogen_idx(water_mol):
    atom_dict = _get_atom_dict(water_mol)
    for k, v in atom_dict.items():
        if v == 1:
            return k


water_mol = _create_water_mol()


def test_jazzy_calibration():
    jazzy_list, _ = calculate_jazzy_and_kallisto_features(water_mol)
    expected = [{'sd': 0.0, 'sa': 1.0, 'num_lp': 2}, {'sd': 1.0, 'sa': 0.0, 'num_lp': 0}, {'sd': 1.0, 'sa': 0.0, 'num_lp': 0}]
    assert jazzy_list == expected


def test_jazzy_encoding():
    jazzy_list, _ = calculate_jazzy_and_kallisto_features(water_mol)
    encoded_list = encode_feature_list('jazzy', jazzy_list)
    oxygen_idx = _get_water_oxygen_idx(water_mol)
    oxygen_vector = encoded_list[oxygen_idx]
    # Check that acceptor strength is not zero
    # by verifying that the first element is equal to 0.
    assert oxygen_vector["sa"][0] == 0
    # Check that donor strength is equal to zero
    # by veryifying that the first element is equal to 1.
    assert oxygen_vector["sd"][0] == 1
    # Check that oxygen has two lone pairs
    assert oxygen_vector["num_lp"][0] == 0
    assert oxygen_vector["num_lp"][2] == 1
    hydrogen_idx = _get_water_hydrogen_idx(water_mol)
    hydrogen_vector = encoded_list[hydrogen_idx]
    # Check that acceptor strength is zero
    # by verifying that the first element is equal to 1.
    assert hydrogen_vector["sa"][0] == 1
    # Check that donor strength is not zero
    # by verifying that the first element is equal to 0.
    assert hydrogen_vector["sd"][0] == 0
    # Check that hydrogen has no lone pairs
    assert hydrogen_vector["num_lp"][0] == 1
    assert hydrogen_vector["num_lp"][1] == 0


def test_kallisto_calibration():
    _, kallisto_list = calculate_jazzy_and_kallisto_features(water_mol)
    expected = [{'eeq': -0.6, 'alp': 5}, {'eeq': 0.3, 'alp': 0}, {'eeq': 0.3, 'alp': 0}]
    assert kallisto_list == expected


def test_kallisto_encoding():
    _, kallisto_list = calculate_jazzy_and_kallisto_features(water_mol)
    encoded_list = encode_feature_list('kallisto', kallisto_list)
    oxygen_idx = _get_water_oxygen_idx(water_mol)
    oxygen_vector = encoded_list[oxygen_idx]
    # Check that eeq is not in either the first or last element.
    assert oxygen_vector["eeq"][0] == 0
    assert oxygen_vector["eeq"][-1] == 0
    # Check that alp is not in the first element.
    assert oxygen_vector["alp"][0] == 0
    hydrogen_idx = _get_water_hydrogen_idx(water_mol)
    hydrogen_vector = encoded_list[hydrogen_idx]
    # Check that eeq is not in either the first or last element.
    assert hydrogen_vector["eeq"][0] == 0
    assert hydrogen_vector["eeq"][-1] == 0
    # Check that alp is in the first element.
    assert hydrogen_vector["alp"][0] == 1
