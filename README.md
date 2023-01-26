# GUTS_MGS_Psych
Analysis scripts GUTS project. Metagenomics shotgun sequencing association with psychiatric patients.

##Software
###Data generation and processing
* (for running .smk file) Snakemake v7.12.0, python v3.1
* (for data cleaning) kneaddata v0.10.0 , Bowtie2 v2.2.5, trimmomatic v0.39-2, human ref for mapping: hg37dec_v0.1
* (taxonomic assingment) motus v3.0.3, python v3.7.15 , SAMtools v1.10, BWA v0.7.17
* (SNV calling and intra-species distance )  motus v3.0.3, python v3.7.15 , SAMtools v1.10, BWA v0.7.17
* (GBM quantifiaction) humann v3.1.0, omixer-rpm v1.1, GBMs database (donwloaded 25-01-2023 from http://raeslab.org/software/gbms.html)
###Data analysis
* (for analysis script) R version v4.0.1
	* General: tidyverse (v1.3.0), vegan (v2.6-4), ape (v5.4), ggrepl (v0.9.1) , MetBrewer (for color) (v0.2.0)
	* Prediction: codacore (0.0.1), tensorflow (2.5.0), pROC (1.17.0.1) 
	* Batch correction: PLSDAbatch (0.2.3)
	* Association: Aldex2 (1.20.0)
	* network formatting: igraph (1.2.5)
* (for bulding the network) SparCC (v0.1.0)
* (for analyzing networks) Netshift webapp

###Datasets
1. GUTS study participants
2. DMP study participants matched to GUTS (age, sex, bmi, socieconomical status) with no dieases
3. (validation 1) Zhu et al. schizophrenia pateitns and healthy controls. ENA accesion: PRJEB29127
4. (validation 2) DMP sutdy participatns with schizophrenia or bipolar disorder





