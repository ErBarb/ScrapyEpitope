# ScrapyEpitope
A pipeline that uses Scrapy and Selenium to access multiple online tools with the aim of predicting and analysing epitopes from a set of proteins.

## How it works
The pipeline uses multiple tools to fulfill its function:

1. `msa.py`: At first it requires a set of proteins as input (swissprot ids), then it uses Multiple Sequence Alignment tools such as MUSCLE and MAFFT to align these sequences. The returned alignment files from these tools are then run through Gblocks (with Selenium) in order to get the set of conserved sequences for each protein within the alignment.

2. `predict.py`: The returned set of sequences from Gblocks is then run through multiple epitope prediction tools, each of which returns a set of predicted epitopes. All the tools, except Ellipro (with Selenium) are accessed by the REST API of IEDB. The tools are:
    * B-Cell Epitope Prediction:
      * Bepipred
      * Emini
      * Kolaskar-Tongaonkar
      * Bepipred2
      * Chou-Fasman
      * Karplus-Schulz
      * Parker
      * Ellipro
    * T-Cell Epitope Prediction:
      * MHC I Binding
      * MHC II Binding
      * MHC I Processing
    * Discontinous Epitope Prediction:
      * Ellipro
      * Discotope(?)

3. `analyse.py`: The sets of epitopes predicted from each of the methods above are then fed into this script, which analyses various parameters for each epitope and stores them in a csv file. The analysis methods and their parameters are:
    *  Protparam - calculates Molecular Weight, Isoelectric Point, Aromaticity, Instability Index, Helix/Turn/Beta Fraction, Molar Ext. Coefficient, Hydropathicity, Flexibility, Charge at pH7
    *  Toxinpred - calculates Toxicity, Hydrophobicity, Steric Hinderance, Side Bulk, Amphipathicity, Hydrophilicity, Net Hydrogen
    *  Algpred - shows whether the epitope is an allergen or not
    *  Vaxijen - shows probability of antigenicity
    *  IEDB Immunogenicity Prediciton Class I and II - Immunogenicity Score
    *  Cluster Analysis - clusters based on sequence identity
    *  Population Coverage (T-Cells only) - calculates the fraction of individuals predicted to respond to a given set of epitopes with known MHC restrictions.
    *  Conservation Analysis - calculates the degree of conservancy of an epitope within a given protein sequence set at different degrees of sequence identity.
    *  Pymol(?) - to visualize the epitopes in the proteins

## How to run the pipeline

  ¯\\\_(ツ)_/¯
