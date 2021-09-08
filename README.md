# Degron Landscape of Fusion Genes in Cancer

This repository contains a pipeline to analyze the impact of gene fusions on protein stability by either degron loss.

## Analysis

As some degrons require a post-translational modification (PTM), you will need to download the PTM data from the PhosphositePlus database (https://www.phosphosite.org/staticDownloads). Once downloaded change the variables PHOS_PATH, ACETYL_PATH and UBIQ_PATH to refer to the location of your download for phosphorylation, acetylation and ubiquitination respectively.

Next you will need to download data from [SNVBox](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3137226/).

```bash
$ cd data
$ wget https://www.dropbox.com/s/gglvpfm06gqzghu/fusion_snvbox_data.zip?dl=1 -O fusion_snvbox_data.zip
$ unzip fusion_snvbox_data.zip
$ cd ..
```

Next, you will need to install the appropriate conda environment (https://docs.conda.io/en/latest/).

```bash
$ conda env create -f environment.yaml
```

Now you need to activate the "fusion" conda environment.

```bash
$ conda activate fusion
```

Running the pipeline is done through the `run.sh` script.

```bash
$ ./run.sh /my/output/dir
```

Note that "/my/output/dir" is whatever directory you choose to output the results as.
