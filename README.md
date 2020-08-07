# MicFunPred:
## More conserved approach to predict functional profiles from 16S sequence data

### Motivation: 
There are multiple tools availbale for the prediction of functional profiles like PICRUSt, Piphillin, Tax4Fun and others. These tools predicts the gene contents of organism without sequenced genomes (at strain level) by using phylogenetic or non-phylogenetic approaches. But, taxonomic identification using one or multiple variable regions of 16S rRNA gene beyond genus level may not be reliable. Secondly, due to database dependency and ever growing reference database size, regular updates are required by these tools. Hence, we proposed a conserved approach to predict gene contents of each genera by considering core genes.

### MicFunPred workflow
MicFunPred relies on ~32,000 genome sequences downloaded from Integrated Microbial Genome database (IMG) representing human, plants, mammals, aquatic and terrestrial ecosystem. 16S rRNA database was constructed using sequences from these genomes and available databases oclustered at 97% identity. MicFunPred is able to predict functional profiles in terms of KEGG Orthology (KO), Enzyme Commission (EC), PFam, TIGRFAM and Cluster of Genes (COG).
![MicFunPred Workflow](workflow.jpeg)

MicFunPred is database/approach independent hence, 16S sequence data processed using QIIME1/2 or DADA2 with any database can be used. MicFunPred follows multiple steps to predict functional profiles:

Input files:

1. Abundance/BIOM table (tab separated)
2. OTUs/ASVs sequences (FASTA format)

#### 1. OTUs/ASVs Sequence mapping
MicFunPred starts by aligning OTUs/ASVs sequences on custom 16S rRNA database using BLAST and assigns genus level taxonomy to each sequence based on user defined percent identity cut-off (-p).

#### 2. Preparation of taxonomy table
The OTU/ASV ids from abundance table are replace by assigned genus and the table is then consolidated and normalized using mean 16S rRNA gene copy numbers of respective genera.

#### 3. Prediction of core genes (Core gene content table)
The core genes for each genera present in abundance table generated above are predicted as per user defined gene coverage cut-off (-c). The core gene is defined as the gene present in x% of the genomes of respective genera. Here, x is gene coverage cut-off and can be adjusted by user (0-1). The value of '0' disables the prediction of core genes.

#### 4. Prediction of metagenomic profiles
The multiplication of abundance table and core gene content table is generated as metagenomic profiles in terms of gene families as stated above.

#### 5. Prediction of pathways
MinPath is used to predict KEGG and MetaCyc pathways with more stringency.

#### 6. (Developmental) Prediction of Community wide Metabolic Potential (CMP)
CMP scores are predicted as suggested in [Noecker C. et al](https://msystems.asm.org/content/1/1/e00013-15) using Predicted Relative Metabolomic Turnover (PRMT) concept.

#### 7. (Developmental) Prediction of phenotypic traits
MicFunPred relies on [BacDive](https://bacdive.dsmz.de/) database for phenotypic information from all bacterial and archeal strains. A seperate BLAST is run to link OTU/ASV sequences with BacDive relitives strictly with 100% identity cut-off.