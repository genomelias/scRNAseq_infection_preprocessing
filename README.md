# scRNAseq_infection_preprocessing


## 1. Data sources

This task is based on publicly available single-cell RNA sequencing data from a study of placental infection. The dataset includes multiple samples unders two different infection conditions (*T. gondii*-infected samples and controls). The single-cell data was originally produces using the Chromium Next GEM Single Cell 3′ v3.1 (Dual index) kit (10x Genomics), sequenced using Illumina Novaseq S4 platform. The sequencing data from 10x Genomics was aligned and quantified using STARsolo v2.7.9a. The subsampled data in Market Exchange Format (MEX) format are stored in `scRNAseq_infection_preprocessing/data/` and are used as the inputs for the workflow.

---

## 2. How to download

The data is stored in ArrayExpress under the identifier E-MTAB-12795. The samples used for this workflow are: Pla_HDBR13007974, Pla_HDBR13007975, Pla_HDBR13798223, and Pla_HDBR13798224. The data is composed by three files per sample: barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz.
### Data can be downloaded as follows:

```bash
wget https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/795/E-MTAB-12795/Files/Pla_HDBR13007974_barcodes.tsv.gz
wget https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/795/E-MTAB-12795/Files/Pla_HDBR13007974_features.tsv.gz
wget https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/795/E-MTAB-12795/Files/Pla_HDBR13007974_matrix.mtx.gz

wget https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/795/E-MTAB-12795/Files/Pla_HDBR13007975_barcodes.tsv.gz
wget https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/795/E-MTAB-12795/Files/Pla_HDBR13007975_features.tsv.gz
wget https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/795/E-MTAB-12795/Files/Pla_HDBR13007975_matrix.mtx.gz

wget https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/795/E-MTAB-12795/Files/Pla_HDBR13798223_barcodes.tsv.gz
wget https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/795/E-MTAB-12795/Files/Pla_HDBR13798223_features.tsv.gz
wget https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/795/E-MTAB-12795/Files/Pla_HDBR13798223_matrix.mtx.gz

wget https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/795/E-MTAB-12795/Files/Pla_HDBR13798224_barcodes.tsv.gz
wget https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/795/E-MTAB-12795/Files/Pla_HDBR13798224_features.tsv.gz
wget https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/795/E-MTAB-12795/Files/Pla_HDBR13798224_matrix.mtx.gz
```


---

## 3. Pre-processing / subsampling

The data was subsampled by only selecting those samples corresponding to placental tissue infected with *Toxoplasma gondii* and their controls at 24 hours post-infections.


---

## 4. How the workflow works
The workflow uses scripts written in Bash and Python. The Bash scripts reformat the downloaded data to make it readable by the Python scripts.
The workflow files is stored in `scRNAseq_infection_preprocessing/workflow/`.

---

### Step 0 – Reformatting data name

**Purpose:** Reformats the downloaded data to make it readable by the Python scripts.
```bash
bash organize_samples.sh
```

---

### Step 1 – Loading data and doublet score calculation

**Purpose:** Loads the data in python and calculated doublet scores.
**Inputs:** MEX files.
**Tools:** `Scanpy`, `Scrublet`
**Outputs:** Scrublet scores.

---

### Step 2 – cell cycle scoring and anndata object creation

**Purpose:** Creates an Anndata object to be used in python. Calculates the cell cycle scores.
**Tools:** `Scanpy`, `anndata`, `numpy`
**Inputs:** MEX files.
**Outputs:** Anndata object.

