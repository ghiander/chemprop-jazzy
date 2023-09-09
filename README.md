# Chemprop-Jazzy
Forked version of the Chemprop framework that allows the calculation of atomistic Jazzy and Kallisto features (https://www.nature.com/articles/s41598-023-30089-x) as described in [Jazzy and Kallisto Features](#jazzy-and-kallisto-features). The functionality currently only works for molecules and not for reactions.

This readme only contains information related to the extension of Chemprop with Jazzy and Kallisto. The official documentation of Chemprop can be found at https://github.com/chemprop/chemprop.

To install the library, refer to the [Installation](#installation) section. Bear in mind that this version of Chemprop can only be installed from source.

Special thanks go Shih-Cheng Li from MIT for his help and support with designing the integration.

## Requirements

For small datasets (~1000 molecules), it is possible to train models within a few minutes on a standard laptop with CPUs only. However, for larger datasets and larger Chemprop models, we recommend using a GPU for significantly faster training.

To use `chemprop` with GPUs, you will need:
 * cuda >= 8.0
 * cuDNN

## Installation

**Note for machines with GPUs:** You may need to manually install a GPU-enabled version of PyTorch by following the instructions [here](https://pytorch.org/get-started/locally/). If you're encountering issues with Chemprop not using a GPU on your system after following the instructions below, check which version of PyTorch you have installed in your environment using `conda list | grep torch` or similar. If the PyTorch line includes `cpu`, please uninstall it using `conda remove pytorch` and reinstall a GPU-enabled version using the instructions at the link above.
̶

### Installing from source (using conda)

1. `git clone https://github.com/ghiander/chemprop-jazzy.git`
2. `cd chemprop-jazzy`
3. `conda env create -f environment.yml`
4. `conda activate jazzprop`
5. `pip install -e .`

### Installing from source (using virtualenv)

1. `git clone https://github.com/ghiander/chemprop-jazzy.git`
2. `cd chemprop-jazzy`
3. `python<version> -m venv jazzprop` (used <version> 3.8 at the time of writing - specifically Python 3.8.5)
4. `source <path_prefix_to_environment>/jazzprop/bin/activate` (e.g. `/home/ghiander/jazzprop/bin/activate`)
5. `python -m pip install flake8 pytest parameterized` (optional - for developers)
6. `python -m pip install -e .` (make sure you are inside `chemprop-jazzy`)

#### Jazzy and Kallisto Features
Jazzy or Kallisto features (https://doi.org/10.1038/s41598-023-30089-x) can be included additionally in the graph convolution by adding the flag `--additional_atom_descriptors` and including the options `jazzy` and/or `kallisto`.

```bash
chemprop_train --data_path tests/data/regression.csv --dataset_type regression --save_dir test_model_checkpoints --quiet --additional_atom_descriptors kallisto jazzy
```

## Examples of usage
Here are some examples of how to use Chemprop-Jazzy

### Training
```bash
# Includes both Kallisto and Jazzy atomic properties
chemprop_train --data_path tests/data/regression.csv --dataset_type regression --save_dir test_model_checkpoints --quiet --additional_atom_descriptors kallisto jazzy

# Includes only Kallisto atomic properties
chemprop_train --data_path tests/data/regression.csv --dataset_type regression --save_dir test_model_checkpoints --quiet --additional_atom_descriptors kallisto

# Includes only Jazzy atomic properties
chemprop_train --data_path tests/data/regression.csv --dataset_type regression --save_dir test_model_checkpoints --quiet --additional_atom_descriptors jazzy

# Includes Jazzy atomic and free hydration energy molecular properties
chemprop_train --data_path tests/data/regression.csv --dataset_type regression --save_dir test_model_checkpoints --quiet --additional_atom_descriptors jazzy --features_generator jazzy_hyd

# Includes Jazzy atomic and hydrogen-bond strength molecular properties
chemprop_train --data_path tests/data/regression.csv --dataset_type regression --save_dir test_model_checkpoints --quiet --additional_atom_descriptors jazzy --features_generator jazzy_hbs
```

### Prediction (via CLI)
```bash
# Chemprop automatically understands whether Jazzy or/and Kallisto were used to train the model
chemprop_predict --test_path tests/data/regression_small.csv --checkpoint_dir test_model_checkpoints --preds_path regression_preds.csv
```

### Prediction (via Python with in-memory model)
```python
import tempfile
import chemprop
import tarfile


def untar_file_to_folder(filepath, output_dir):
    with tarfile.open(filepath) as tf:
        tf.extractall(output_dir)


# Extract the model from an archive
model_tar = "tests/data/regression_model_jazzy_kallisto.tar.gz"
with tempfile.TemporaryDirectory() as tmp_folder:
    untar_file_to_folder(model_tar, tmp_folder)

    # Configure arguments
    args_list = ['--test_path', 'tests/data/regression_small.csv',
                 '--preds_path', '/dev/null',
                 '--checkpoint_dir', tmp_folder]
    args = chemprop.args.PredictArgs().parse_args(args_list)

    # Load model in memory
    model_objects = chemprop.train.load_model(args=args)

    # Predict values for a list of SMILES
    smiles_list = [['CCC'], ['CCCC'], ['OCC']]
    print(chemprop.train.make_predictions(args=args,
                                          smiles=smiles_list,
                                          model_objects=model_objects))
```

```
>>> ...
[[-3.6102219301449665], [-3.5861583065877727], [-3.361401796920642]]
```

### For developers
To run the test suite use `pytest -v`
