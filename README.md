# acidify-ME
A ME-model framework to simulate the response of E. coli under acid stress

### Requirements

* python == 2.7
* pandas >= 0.21.0
* numpy >= 1.13.3
* cobrame >= 0.0.8
* cobrapy == 0.5.11
* biopython >= 1.60

For more information on solving the ME-model, please refer to [cobrame package](https://github.com/SBRG/cobrame). The related [paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007525) has been published at PLOS Computational Biology.

## Overview of files in /data
* /proteins, a directory that contains the information about the periplasmic proteins of *E. coli* as described by the ME-model, including protein structure (prot.pdb) and the calculated charges of the folded protein and its side chains under different pH values (sum_crg.out).
* b3509.fasta, sequence of the *hdeB* gene.
* charged_aa_side_chain_pKa.csv, table of individual amino acid side chain pKa values, used to calculate the charge of unfolded protein (peptide) at different pH values.
* lipid_fatty_acid_composition_under_acid.csv, relative fraction of membrane lipids with different fatty acid tails in terms of non-adapted and acid-adapted profiles.
* peptide_radius_of_gyration.csv, table of empirical measurements on the radius of gyration of peptides.
* periplasm_protein_fold_rate.csv, fold rates of periplasmic proteins, more information on how protein fold rate is calculated can be found in [ssbio package](https://github.com/SBRG/cobrame).
