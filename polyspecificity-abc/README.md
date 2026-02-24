# ABC Polyspecificity

This folder contains code and example data to run SVM (support vector machine) predictions for our publication:

**A comprehensive protein structures- and machine learning-based workflow to identify polyspecific binding sites for drug development**

by Gergely Gyimesi, Shouyu Wang, Hauke Busch, Sven Marcel Stefan.

## Installation

1. Create a `conda` environment based on `environment.yaml`, which should also install `gemmi` and `reduce`.

```bash
conda env create -f environment.yaml -p ./pyenv
```

2. Download and install MSMS into `/opt/msms` from https://ccsb.scripps.edu/msms/

3. Download and install MGLTools from https://ccsb.scripps.edu/mgltools/downloads/. The expected path is `/opt/mgltools_x86_64Linux2_1.5.7`. First untar the archive and move it into the destination path (you might have to use `sudo`). Then, `cd ` into that directory and use the `./install.sh` script in MGLTools without parameters to install into the current directory. This way, it will also install the `python2.7` binary into the `bin/` folder.

4. You will need GNU make to run the calculations. Make sure that GNU make is installed, for example, on a Debian/Ubuntu machine:

```bash
sudo apt install make
```

The installation and usage have been tested and should work on Ubuntu 24.04 LTS.

## Usage

### PDB repository

Structures for analysis should be placed into a PDB repository, as an example, the `test-pdb` directory is provided with some example structures. The structures are organized in a similar way like in the PDB (*i.e.*, in subfolders according to the second and third characters of the PDB ID). Here, structures in mmCIF format are expected in the `test-pdb/mmcif` folder, if you want to use biological assemblies, place them into the `test-pdb/biounit-mmcif` folder according to the examples.

For each structure, you need to generate mmJSON files:
```bash
cd test-pdb
make all-mmjsons all-biounit-mmjsons
cd ..
```

### Generating complexes

Proximity information needs to be generated first for all structures:
```bash
bash gen-all-proximities.sh
```

This command will generate a `proximity.5.tsv` storing information about which molecules are close to which other molecules within the structures.

At this point, you need to create or edit the `complexes.tsv` file, and include basic information about the ligands that you would like to study. The file is a text file containing tab-separated values, be sure to preserve the tab characters when editing. Ligands are identified by the PDB ID of the structure (`pdb_id` column), their 3-letter or 5-letter PDB code (`mon_id` column), and the asymmetric ID within the mmCIF file (`asym_id` column, also referred to as `label_asym_id` in the `_atom_site` records of mmCIF files). See the existing file for example values.

Once the `complexes.tsv` file is final, generate the complex structures and other derived files for each ligand instance:
```bash
python gen-complex-strucs.py complexes.tsv proximity.5.tsv
make -C complex-strucs all
```

### Surface area calculation

The next step is to calculate surface areas of the protein around each ligand instance, as well as the surface area of the ligands themselves.

```bash
make -C complex-strucs/surface-area all
```

### Running SVM predictions

Once the surface areas are calculated, you are ready to calculate normalized surface areas required for the SVM to process. These will be stored into the tabular file `surfareas.tsv`:
```bash
python get-all-surfareas.py >surfareas.tsv
```

Finally, to perform predictions, run the following command, which will apply the best-performing SVM as described in the publication:
```bash
python svm-predict.py surfareas.tsv svm3.best.pkl >svm-predictions.tsv
```

The resulting `svm-predictions.tsv` is a tab-separated file that can be loaded into spreadsheet programs, and should contain a `svc_pred` column with predictions (currently either "on-target" or "off-target").

