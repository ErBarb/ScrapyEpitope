# ScrapyEpitope
A pipeline that uses Scrapy and Selenium to access multiple online tools with the aim of predicting and analysing epitopes from a set of proteins.

## How it works
The pipeline uses multiple tools to fulfill its function:

1. msa.py: In the beginning it requires a set of proteins as input (swissprot ids), whereby it uses Multiple Sequence Alignment tools such as MUSCLE and MAFFT to align these sequences. The returned alignment files from these tools are then run through Gblocks (Selenium) in order to get the set of conserved sequences for each protein within the alignment.

2. predict.py: The returned set of sequences from Gblocks is then run through multiple epitope prediction tools, each of which returns a set of predicted epitopes. The tools are:
  * B-Cell Epitope Prediction:
    * Bepipred
    * Emini
    * Kolaskar-Tongaonkar
    * Bepipred2
    * Chou-Fasman
    * Karplus-Schulz
    * Parker
  * T-Cell Epitope Prediction:
    * IEDB MHC Class I
    * IEDB MHC Class II
    * NetCTLpan
  * Discontinous Epitope Prediction:
    * Ellipro
    * Discotope

3. analyse.py: The sets of epitopes predicted from each of the methods above are then fed into this script, which analyses various parameters for each epitope and stores them in a csv file. The analysis methods are:
  *  Protparam
  *  Toxinpred
  *  Algpred
  *  Vaxijen
  *  IEDB Immunogenicity Prediciton Class I and II - only for T-Cells
  *  Cluster Analysis
  *  Population Coverage - only for T-Cells
  *  Conservation Analysis
  *  Pymol - visualization

## How to run the pipeline

  ¯\\\_(ツ)_/¯
