# phage_typing
https://www.frontiersin.org/articles/10.3389/fmicb.2018.00836/full

## Description

Prophage typing tool from genome sequences. Main steps are:

* Basic assembly of Illumina paired-end reads with SPAdes
* Detection of the prophage with PHASTER
* Clustering of prophage sequences with CD-HIT-EST
* Prophage diversity anlysis with QIIME

## Requirements

* QIIME. Best is to create a virtual environment, as recommended by the authors.
* CD-HIT-EST
* BBTools
* SPAdes
* Java
* Python 3
* GNU parallel
* Pigz
* Biom

## Usage
It is recommended to use the "phage_diversity_HT.sh" script instead of "phage_diversity.sh", where the "HT" stands for "High Throughput". The "HT" version do all the read pre-processing and assembling in a single fonction, and the temporary files are deleted when no longer needed, which allow to process a high number of samples with limited hard drive space. The trade off is that version has less QC included. The "HT" version is also the most up-to-date one.

In its current implementation, it is better to run the script by executing sections at the time, i.e. copy pasting from script to terminal. 
