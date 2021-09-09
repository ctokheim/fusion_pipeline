# Degron Landscape of Fusion Genes in Cancer

This repository contains a pipeline to analyze the impact of gene fusions on protein stability by either degron loss.

## Download miscellaneous data

As some degrons require a post-translational modification (PTM), you will need to download the PTM data from the PhosphositePlus database (https://www.phosphosite.org/staticDownloads). Once downloaded change the variables PHOS_PATH, ACETYL_PATH and UBIQ_PATH to refer to the location of your download for phosphorylation, acetylation and ubiquitination respectively.

Next you will need to download data from [SNVBox](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3137226/).

```bash
$ cd data
$ wget https://www.dropbox.com/s/gglvpfm06gqzghu/fusion_snvbox_data.zip?dl=1 -O fusion_snvbox_data.zip
$ unzip fusion_snvbox_data.zip
$ cd ..
```

## Install environment

This code only works on *linux* and was tested on CentOS release 6.10. No non-standard hardware are needed. The code has only been tested on CentOS 6.10 and software version numbers listed in the environment.yaml config.

Next, you will need to install the appropriate conda environment (https://docs.conda.io/en/latest/). 

```bash
$ conda env create -f environment.yaml
```

Note that software dependencies are expressed in the environment.yaml config file. It may take ~30min to 1 hour to install dependencies. Now you need to activate the "fusion" conda environment.

```bash
$ conda activate fusion
```

## Analysis

Running the pipeline is done through the `run.sh` script. Note that by default it uses the fusion calls that were analyzed in the paper.

```bash
$ ./run.sh output
```

Note that you can change "output" to whatever directory you choose to output the results as.

You will need to stop midway through the shell script (line #47), and then use the "motif/wt_motif_hits_canonical.snvbox_input.txt" file as input to [SNVBox](https://wiki.chasmsoftware.org/index.php/SNVBox_Tutorial). Once finished, copy the file to "degron_pred/wt_motif_hits_canonical.snvbox_annot.txt" location within your output directory and run the remaining portions of the shell script.

Running all of the commands may take ~10-20 hours to run on a normal computer without parallelization.

Ultimately, you should get a similar result as [this](https://www.dropbox.com/s/8uek4y2lflowhrq/example_result.tar.gz?dl=1).
