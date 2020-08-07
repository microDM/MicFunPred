#!/usr/bin/python

'''
Created on 22-Feb-2019

@author: Dattatray Mongad,
                    Junior Research Fellow (PhD Student),
                    National Centre fo Microbial Resource,
                    National Centre for Cell Science,
                    Pune, Maharshtra, India.
'''
import argparse
from micfunpreDefinitions import micfunpreDefinitions as fp
import pandas as pd
import warnings
import os
import shutil
import time
import subprocess

# make options
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--otu_table", metavar="PATH",type=str,
                  help="[Required] Tab delimited OTU table")
parser.add_argument("-r", "--repset_seq",metavar="PATH",type=str,
                  help="[Required] Multi-fasta file of OTU/ASV sequences")
parser.add_argument("-p", "--perc_ident",type=float,
                  help="[Optional] Percent identity cut-off to assign genus. (default: 97)", default=98.3)
parser.add_argument("-b","--blastout",type=str,help="Blast output of ASV/OTU sequences with any database in output format 6",
                    default=None,metavar="PATH")
parser.add_argument("-c", "--coreper",type=str, default=0.5,
                  help="[Optional] Percentage of organism in a genus which should have gene to define it as core. Value ranges "
                       "from 0 to 1 (default: 0.5)")
parser.add_argument("-o", "--output", type=str,metavar="PATH",
                  help="[Optional] Output directory (default: funpred_out)",
                  default="funpred_out")
parser.add_argument("-l",'--tax_level',default='genus',type=str,
                  help='[Optional] (genus OR species). Taxonomy level to which analysis to be done. If '
                  'species is given, OTUs/ASVs which show 100%% identity to database sequence will be given '
                  'species level taxonomy')
parser.add_argument("--predict_trait",default=False,
                  help='[Optional] Predict phenotypic metadata also. This step may take some time.(default: 0)',
                  action='store_true')
parser.add_argument("-t", "--threads",type=str,metavar="INT",
                  help="(Optional) number of threads to be used. (default: 1)", default="1")
parser.add_argument("-v", "--verbose", action="store_true", default=False,
                  help="Print message of each step to stdout.")

# if required option does not provided exit
options = parser.parse_args()
if (options.otu_table and options.repset_seq) is None:
    exit("Please provide all required inputs or type python3.6 runPrediction.py -h")
else:
    in_otuTab = options.otu_table
    in_repSet = options.repset_seq
    in_out = options.output
    in_numCores = options.threads
    in_percIdentCutOff = options.perc_ident
    in_coreperc = float(options.coreper)
    in_verbose = options.verbose
    in_phenotype = options.predict_trait
    in_blastout = options.blastout

###########################################        MAIN         ######################33

#data path
funpredPath = fp.baseDir()

cwd = os.getcwd()
os.chdir(cwd)

dbPath = funpredPath + "/data/db/"
dataPath = funpredPath + "/data/"
otherPath = funpredPath + "/data/other/"

# handle if directory already exists
if (os.path.exists(in_out)):
    flag = input("Directory " + in_out + " already exists. Enter \"1\" to overwrite or \"0\" to exit. \n")
    if (flag == "1"):
        shutil.rmtree(in_out)
        os.mkdir(in_out)
        cwd = cwd + "/" + in_out
    else:
        exit("Exiting...")
else:
    os.mkdir(in_out)
    cwd = cwd + "/" + in_out

# calculate time
start_time = time.time()

if (in_verbose):
    print("Running BLAST with " + str(in_percIdentCutOff) + " cut-off.")
if (in_blastout == None):
    fp.blast(os.path.join(os.getcwd(),in_repSet), dbPath + "blastDB/16S_all", in_numCores, cwd)
    tax_dict, abundTable = fp.selectBlastHits_assignGenus_subsetOtuTable(cwd + "/out.blast", in_otuTab,
                                                                     in_percIdentCutOff,in_level)
elif(in_blastout != None):
    #filter with pident cut-off
    df = pd.read_csv(in_blastout,sep='\t',header=None)
    in_blastout = cwd+'/'+in_blastout+'_'+str(in_percIdentCutOff)
    df.loc[df[2]>=float(in_percIdentCutOff)].to_csv(in_blastout,sep='\t',header=None,index=None)
    tax_dict, abundTable = fp.selectBlastHits_assignGenus_subsetOtuTable(in_blastout, in_otuTab,
                                                                     in_percIdentCutOff,in_level)

abundTable.to_csv(cwd + "/tax_abund.table", sep="\t")

if (in_verbose):
    print("Predicting 16S rRNA copy numbers")
# read 16S rRNA copy number table
copyNumberTable_16S = pd.concat([pd.read_csv(x,sep='\t',index_col=0) for x in [os.path.join(dbPath,x) for x in os.listdir(dbPath) if('16S_cp' in x)]],sort=False)
copyNumberTable_16S_consolidated = fp.makeTable16S(copyNumberTable_16S, "mean", list(abundTable.index))
copyNumberTable_16S_consolidated.to_csv(cwd + "/predicted_16S_copy_numbers.txt", sep="\t", header=None)

# divide abundance table by 16S rRNA copy number
for index in list(copyNumberTable_16S_consolidated.index):
    abundTable.loc[index] = abundTable.loc[index] / float(copyNumberTable_16S_consolidated.loc[index])
abundTable = abundTable.round(2).fillna(0)
abundTable.to_csv(cwd + "/tax_abund_normalized.table", sep="\t")

# KEGG prediction
# 1. Predict gene copy numbers
if (in_verbose):
    print("Predicting KO copy numbers")
copyNumberTable_gene = pd.read_parquet(os.path.join(dbPath,'ko.parq'))
copyNumberTable_gene_consolidated = fp.makeKOTable(copyNumberTable_gene, abundTable, in_coreperc).fillna(0)
del copyNumberTable_gene # remove large table to save memory
copyNumberTable_gene_consolidated.to_csv(cwd + "/predicted_KO.txt", sep="\t")
# 2. Multiplication
if (in_verbose):
    print("Predicting KO metagenome")
final_df = abundTable.transpose().dot(copyNumberTable_gene_consolidated).transpose()
final_df = final_df.loc[final_df.sum(axis=1) != 0]
os.mkdir(cwd + "/KO_metagenome")
final_df.to_csv(cwd + "/KO_metagenome/KO_metagenome.txt", sep="\t")
del copyNumberTable_gene_consolidated
# 3. Add description
if (in_verbose):
    print("Anotating predicted metagenome")
final_df = fp.addAnnotations(final_df, otherPath + "ko00001.txt")
final_df.to_csv(cwd + "/KO_metagenome/KO_metagenome_with_description.txt", sep="\t")
# 4. MinPath
if (in_verbose):
    print("Running MinPath (KO)")
final_df = fp.runMinPath(final_df, funpredPath, cwd + "/KO_metagenome", "kegg")
# 5. Groupby levels
fp.summarizeByFun(final_df, "A").to_csv(cwd + "/KO_metagenome/summarized_by_A.txt", sep="\t")
fp.summarizeByFun(final_df, "B").to_csv(cwd + "/KO_metagenome/summarized_by_B.txt", sep="\t")
fp.summarizeByFun(final_df, "C").to_csv(cwd + "/KO_metagenome/summarized_by_C.txt", sep="\t")
fp.summarizeByFun(final_df, "Pathway_Module").to_csv(cwd + "/KO_metagenome/summarized_by_Pathway_Module.txt",
                                                             sep="\t")
del final_df

# EC prediction
# 1. Predict gene copy numbers
if (in_verbose):
    print("Predicting EC copy numbers")
copyNumberTable_gene = pd.read_parquet(os.path.join(dbPath,'ec.parq'))
copyNumberTable_gene_consolidated = fp.makeKOTable(copyNumberTable_gene, abundTable, in_coreperc).fillna(0)
del copyNumberTable_gene # remove large table to save memory
copyNumberTable_gene_consolidated.to_csv(cwd + "/predicted_EC.txt", sep="\t")
# 2. Multiplication
if (in_verbose):
    print("Predicting EC metagenome")
final_df = abundTable.transpose().dot(copyNumberTable_gene_consolidated).transpose()
final_df = final_df.loc[final_df.sum(axis=1) != 0]
os.mkdir(cwd + "/MetaCyc_metagenome")
final_df.to_csv(cwd + "/MetaCyc_metagenome/EC_metagenome.txt", sep="\t")
del copyNumberTable_gene_consolidated
# 3. MinPath
if (in_verbose):
    print("Running MinPath (RXN)")
final_df = fp.ec2RXN(final_df, otherPath + "ec2rxn_new")
final_df.to_csv(cwd + "/MetaCyc_metagenome/RXN_metagenome.txt", sep="\t")
os.mkdir(cwd + "/MetaCyc_metagenome" + "/minPath_files")
final_df = fp.runMinPath(final_df, funpredPath, cwd + "/MetaCyc_metagenome/minPath_files",
                                          "metacyc")
final_df.to_csv(cwd + "/MetaCyc_metagenome/PathwayAbundance.table", sep="\t")
# 4. Groupby levels
final_df = fp.addMetaCycPathwayName(final_df, otherPath + "path_to_Name.txt")
final_df.to_csv(cwd + "/MetaCyc_metagenome/PathwayAbundance_with_names.table", sep="\t")
fp.summarizeByFun(final_df, "Type").to_csv(
    cwd + "/MetaCyc_metagenome/Pathway_summarize_by_Types.table", sep="\t")
del final_df

# Pfam prediction
# 1. Predict gene copy numbers
if (in_verbose):
    print("Predicting Pfam copy numbers")
copyNumberTable_gene = pd.read_parquet(os.path.join(dbPath,'pfam.parq'))    
copyNumberTable_gene_consolidated = fp.makeKOTable(copyNumberTable_gene, abundTable, in_coreperc).fillna(0)
del copyNumberTable_gene
copyNumberTable_gene_consolidated.to_csv(cwd + "/predicted_Pfam.txt", sep="\t")
# 2. Multiplication
final_df = abundTable.transpose().dot(copyNumberTable_gene_consolidated).transpose().round()
final_df = final_df.loc[final_df.sum(axis=1) != 0]
os.mkdir(cwd + "/Pfam_metagenome")
final_df.to_csv(cwd + "/Pfam_metagenome/Pfam_metagenome.txt", sep="\t")
del copyNumberTable_gene_consolidated
# 3. Add description
final_df = fp.addAnnotations(final_df, otherPath + "pfam_mapping.txt")
final_df.to_csv(cwd + "/Pfam_metagenome/Pfam_metagenome_with_description.txt", sep="\t")
del final_df

# COG prediction
# 1. Predict gene copy numbers
if (in_verbose):
    print("Predicting COG copy numbers")
copyNumberTable_gene = pd.read_parquet(os.path.join(dbPath,'cog.parq'))
copyNumberTable_gene_consolidated = fp.makeKOTable(copyNumberTable_gene, abundTable, in_coreperc).fillna(0)
del copyNumberTable_gene
copyNumberTable_gene_consolidated.to_csv(cwd + "/predicted_COG.txt", sep="\t")
# 2. Multiplication
final_df = abundTable.transpose().dot(copyNumberTable_gene_consolidated).transpose().round()
final_df = final_df.loc[final_df.sum(axis=1) != 0]
os.mkdir(cwd + "/COG_metagenome")
final_df.to_csv(cwd + "/COG_metagenome/COG_metagenome.txt", sep="\t")
del copyNumberTable_gene_consolidated
# 3. Add description
final_df = fp.addAnnotations(final_df, otherPath + "cognames2003-2014.tab")
final_df.to_csv(cwd + "/COG_metagenome/COG_metagenome_with_description.txt", sep="\t")
del final_df

# TIGRFAM prediction
# 1. Predict gene copy numbers
if (in_verbose):
    print("Predicting TIGRFAM copy numbers")
copyNumberTable_gene = pd.read_parquet(os.path.join(dbPath,'tigrfam.parq'))
copyNumberTable_gene_consolidated = fp.makeKOTable(copyNumberTable_gene, abundTable, in_coreperc).fillna(0)
del copyNumberTable_gene
copyNumberTable_gene_consolidated.to_csv(cwd + "/predicted_TIGRFAM.txt", sep="\t")
# 2. Multiplication
final_df = abundTable.transpose().dot(copyNumberTable_gene_consolidated).transpose().round()
final_df = final_df.loc[final_df.sum(axis=1) != 0]
os.mkdir(cwd + "/TIGRFAM_metagenome")
final_df.to_csv(cwd + "/TIGRFAM_metagenome/TIGRFAM_metagenome.txt", sep="\t")
del copyNumberTable_gene_consolidated
# 3. Add description
final_df = fp.addAnnotations(final_df, otherPath + "TIGRFAMs_9.0_INFO.txt")
final_df.to_csv(cwd + "/TIGRFAM_metagenome/TIGRFAM_metagenome_with_description.txt", sep="\t")
del final_df

timeTaken = round(time.time() - start_time, 2)
print("Completed pipline in " + str(timeTaken) + " seconds.")



