import logging
from chemprop.features.encoding import bin_fraction
from chemprop.features.encoding import bin_value
from jazzy.core import kallisto_molecule_from_rdkit_molecule
from jazzy.core import get_covalent_atom_idxs
from jazzy.core import get_charges_from_kallisto_molecule
from jazzy.core import calculate_polar_strength_map
from jazzy.utils import JazzyError
from rdkit import Chem
from rdkit.Chem import AllChem


logging.basicConfig(
    format="Jazzy %(levelname)s: [%(asctime)s] %(message)s",
    level=logging.WARNING,
    datefmt="%H:%M:%S",
)


def calculate_jazzy_and_kallisto_features(mh):
    """
    Calculates lists of Jazzy and Kallisto dictionaries.
    One dictionary per atom is created.
    """
    mh = _embed_and_minimise_rdkit_molecule(mh)
    kallisto_mol = kallisto_molecule_from_rdkit_molecule(mh)
    atoms_and_nbrs = get_covalent_atom_idxs(mh)
    kallisto_charges = get_charges_from_kallisto_molecule(kallisto_mol,
                                                          charge=0)
    pol_map = calculate_polar_strength_map(mh, kallisto_mol,
                                           atoms_and_nbrs,
                                           kallisto_charges)
    jazzy_list, kallisto_list = _process_features(pol_map)
    return jazzy_list, kallisto_list


def _embed_and_minimise_rdkit_molecule(mh):
    """Generate coordinates and minimise energy using MMFF94 for a molecule."""
    mh = _embed_rdkit_molecule(mh)
    AllChem.MMFFOptimizeMolecule(mh)
    return mh


def _embed_rdkit_molecule(mh):
    """Generates 3D coordinates for a molecule."""
    emb_code = AllChem.EmbedMolecule(mh, randomSeed=11)
    if emb_code == -1:
        logging.error("The RDKit embedding has failed "
                      f"for the molecule: {Chem.MolToSmiles(mh)}")
        raise JazzyError("Could not embed the compound.")
    return mh


def _process_features(pol_map):
    """Extracts atomic features into two lists of atom dicts."""
    jazzy_list = list()
    kallisto_list = list()
    for feat_dict in pol_map.values():
        jazzy_list.append(__get_jazzy_features(feat_dict))
        kallisto_list.append(__get_kallisto_features(feat_dict))
    return jazzy_list, kallisto_list


def __get_jazzy_features(feat_dict):
    """
    Selects and condenses Jazzy features.
    Rounding is needed for binning.
    """
    return {
            "sd": bin_fraction(feat_dict["sdc"] + feat_dict["sdx"], 0.1),
            "sa": bin_fraction(feat_dict["sa"], 0.1),
            "num_lp": feat_dict["num_lp"]
        }


def __get_kallisto_features(feat_dict):
    """
    Selects and rounds Kallisto features.
    Rounding is needed for binning. Alp is
    rounded to the nearest 5.
    """
    return {
            "eeq": bin_fraction(feat_dict["eeq"], 0.1),
            "alp": bin_value(feat_dict["alp"], 5)
        }
