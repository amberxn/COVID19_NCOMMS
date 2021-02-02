# COVID19_NCOMMS

This repository includes custom codes to perform metabolomics data normalization, process cytokine data, perform regression analysis, and cluster cytokines and metabolits by Mfuzz in **Host metabolite-cytokine correlation landscape in SARS-CoV-2 infection**.

#### Normalization_QCmad-TIC.R
To perform QCmad-TIC normalization of targeted and untargeted metabolomics NA-replaced raw peak area data.
> Metabolomics raw peak area data are provided as Supplementary Data 2 with this paper.\
> Normalized metabolomics data are provided as Supplementary Data 3 with this paper.\
> Source data are provided in Source Data file with this paper.
#### MetabCyto_process.R
To analyse significanly altered metabolites and cytokines, perform KEGG metabolic pathway enrichemnt analysis. The related figures are Fig. 1b, e and Supplementary Fig. 2e-f; 3a.
> Source data are provided in Source Data file with this paper.
#### CytoMetab_regression.R
To perform regression analysis between cytokines and metabolites in severe patients, mild patients, and healthy controls, separately, and to perform subsequent KEGG metabolic pathway enrichment analysis and unsupervised clustering by igraph. The related figures are Fig. 2 and Supplementary Fig. 4.
> Source data are provided in Source Data file with this paper.
#### Longitudinal_regression.R
To cluster cytokines and mentaolites by Mfuzz, to perform regression analysis between cytokines and metabolites in follow-up patients, and to perform subsequent KEGG metabolic pathway enrichment analysis. The related figures are Fig. 3 and Supplementary Fig. 5; 6.
> Source data are provided in Source Data file with this paper.

#### Citation
Host metabolite-cytokine correlation landscape in SARS-CoV-2 infection. *Nat Commun* (2021). 

Thanks for your interest :)
