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
                  help="[Optional] Percent identity cut-off to assign genus. (default: 98.3)", default=98.3)
parser.add_argument("-b","--blastout",type=str,help="Blast output of ASV/OTU sequences with any database in output format 6",
                    default=None,metavar="PATH")
parser.add_argument("-c", "--coreper",type=str, default=1,
                  help="[Optional] Percentage of organism in a genus which should have gene to define it as core. Value ranges "
                       "from 0 to 1 (default: 1)")
parser.add_argument("-o", "--output", type=str,metavar="PATH",
                  help="[Optional] Output directory (default: funpred_out)",
                  default="funpred_out")
parser.add_argument("-l",'--tax_level',default='genus',type=str,
                  help='[Optional] (genus OR species). Taxonomy level to which analysis to be done. If '
                  'species is given, OTUs/ASVs which show 100%% identity to database sequence will be given '
                  'species level taxonomy')
choice_db = ("all","human","plants","aquatic","mammals","terrestrial")
parser.add_argument("-d",'--db',default="all",help="From which environment database genomes should be used",
                    choices=choice_db)
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
    in_level = options.tax_level
    in_blastout = options.blastout
    in_db = options.db

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
    fp.blast(in_repSet, dbPath + "blastDB/16S_all", in_numCores, cwd)
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
#copyNumberTable_16S,copyNumberTable_ko,copyNumberTable_ec = "","",""

if (in_verbose):
    print("Predicting 16S rRNA copy numbers")
# read 16S rRNA copy number table
if(in_db == "all"):
    copyNumberTable_16S = pd.concat([pd.read_csv(x,sep='\t',index_col=0) for x in [os.path.join(dbPath,x) for x in os.listdir(dbPath) if('16S_cp' in x)]],sort=False)    
else:
    copyNumberTable_16S = pd.read_csv(dbPath + in_db + "_16S_cp.txt", sep="\t", index_col=0)
copyNumberTable_16S_consolidated = fp.makeTable16S(copyNumberTable_16S, "mean", list(abundTable.index))
copyNumberTable_16S_consolidated.to_csv(cwd + "/predicted_16S_copy_numbers.txt", sep="\t", header=None)

# divide abundance table by 16S rRNA copy number
for index in list(copyNumberTable_16S_consolidated.index):
    abundTable.loc[index] = abundTable.loc[index] / float(copyNumberTable_16S_consolidated.loc[index])
abundTable = abundTable.round(2).fillna(0)
abundTable.to_csv(cwd + "/tax_abund_normalized.table", sep="\t")

# make ko table for given taxonomy
if (in_verbose):
    print("Predicting KO copy numbers")
if(in_db == "all"):
    copyNumberTable_ko = pd.concat([pd.read_csv(x,sep='\t',index_col=0, compression='gzip') for x in [os.path.join(dbPath,x) for x in os.listdir(dbPath) if('ko.tsv.gz' in x)]],sort=False)  
else:
    copyNumberTable_ko = pd.read_csv(dbPath + in_db + "_ko.tsv.gz", sep="\t", index_col=0, compression='gzip')
df_cutoff = pd.read_csv(funpredPath+'/data/GC_KO.txt',sep='\t',index_col=0)
copyNumberTable_ko_consolidated = fp.makeKOTable(copyNumberTable_ko, abundTable, in_coreperc,df_cutoff).fillna(0)
copyNumberTable_ko = None  # remove large table to save memory
copyNumberTable_ko_consolidated.to_csv(cwd + "/predicted_KO.txt", sep="\t")

"""
if (in_verbose):
    print("Predicting EC copy numbers")
if(in_db == "all"):
    copyNumberTable_ec = pd.concat([pd.read_csv(x,sep='\t',index_col=0, compression='gzip') for x in [os.path.join(dbPath,x) for x in os.listdir(dbPath) if('ec.tsv.gz' in x)]],sort=False)    
else:
    copyNumberTable_ec = pd.read_csv(dbPath + in_db + "_ec.tsv.gz", sep="\t", index_col=0, compression="gzip")
# make EC table for given taxonomy
copyNumberTable_ec_consolidated = fp.makeKOTable(copyNumberTable_ec, abundTable, in_coreperc).fillna(0)
copyNumberTable_ec = None
copyNumberTable_ec_consolidated.to_csv(cwd + "/predicted_EC.txt", sep="\t")

if (in_verbose):
    print("Predicting Pfam copy numbers")
if(in_db == "all"):
    copyNumberTable_pfam = pd.concat([pd.read_csv(x,sep='\t',index_col=0, compression='gzip') for x in [os.path.join(dbPath,x) for x in os.listdir(dbPath) if('pfam.tsv.gz' in x)]],sort=False)    
else:
    copyNumberTable_pfam = pd.read_csv(dbPath + in_db + "_pfam.tsv.gz", sep="\t", index_col=0, compression="gzip")
# make EC table for given taxonomy
copyNumberTable_pfam_consolidated = fp.makeKOTable(copyNumberTable_pfam, abundTable, in_coreperc).fillna(0)
copyNumberTable_pfam = None
copyNumberTable_pfam_consolidated.to_csv(cwd + "/predicted_Pfam.txt", sep="\t")

if (in_verbose):
    print("Predicting COG copy numbers")
if(in_db == "all"):
    copyNumberTable_cog = pd.concat([pd.read_csv(x,sep='\t',index_col=0, compression='gzip') for x in [os.path.join(dbPath,x) for x in os.listdir(dbPath) if('cog.tsv.gz' in x)]],sort=False)    
else:
    copyNumberTable_cog = pd.read_csv(dbPath + in_db + "_cog.tsv.gz", sep="\t", index_col=0, compression="gzip")
# make EC table for given taxonomy
copyNumberTable_cog_consolidated = fp.makeKOTable(copyNumberTable_cog, abundTable, in_coreperc).fillna(0)
copyNumberTable_cog = None
copyNumberTable_cog_consolidated.to_csv(cwd + "/predicted_COG.txt", sep="\t")

if (in_verbose):
    print("Predicting TIGRFAM copy numbers")
if(in_db == "all"):
    copyNumberTable_tigrfam = pd.concat([pd.read_csv(x,sep='\t',index_col=0, compression='gzip') for x in [os.path.join(dbPath,x) for x in os.listdir(dbPath) if('tigrfam.tsv.gz' in x)]],sort=False)    
else:
    copyNumberTable_tigrfam = pd.read_csv(dbPath + in_db + "_tigrfam.tsv.gz", sep="\t", index_col=0, compression="gzip")
# make EC table for given taxonomy
copyNumberTable_tigrfam_consolidated = fp.makeKOTable(copyNumberTable_tigrfam, abundTable, in_coreperc).fillna(0)
copyNumberTable_tigrfam = None
copyNumberTable_tigrfam_consolidated.to_csv(cwd + "/predicted_TIGRFAM.txt", sep="\t")

"""

# multiplication
# KO
if (in_verbose):
    print("Predicting KO, EC, Pfam, COG and TIGRFAM metagenome")
final_df = abundTable.transpose().dot(copyNumberTable_ko_consolidated).transpose()
final_df = final_df.loc[final_df.sum(axis=1) != 0]
os.mkdir(cwd + "/KO_metagenome")
final_df.to_csv(cwd + "/KO_metagenome/KO_metagenome.txt", sep="\t")

"""
# EC
final_df_ec = abundTable.transpose().dot(copyNumberTable_ec_consolidated).transpose().round()
final_df_ec = final_df_ec.loc[final_df_ec.sum(axis=1) != 0]
os.mkdir(cwd + "/MetaCyc_metagenome")
final_df_ec.to_csv(cwd + "/MetaCyc_metagenome/EC_metagenome.txt", sep="\t")
# Pfam
final_df_pfam = abundTable.transpose().dot(copyNumberTable_pfam_consolidated).transpose().round()
final_df_pfam = final_df_pfam.loc[final_df_pfam.sum(axis=1) != 0]
os.mkdir(cwd + "/Pfam_metagenome")
final_df_pfam.to_csv(cwd + "/Pfam_metagenome/Pfam_metagenome.txt", sep="\t")
# COG
final_df_cog = abundTable.transpose().dot(copyNumberTable_cog_consolidated).transpose().round()
final_df_cog = final_df_cog.loc[final_df_cog.sum(axis=1) != 0]
os.mkdir(cwd + "/COG_metagenome")
final_df_cog.to_csv(cwd + "/COG_metagenome/COG_metagenome.txt", sep="\t")
# TIGRFAM
final_df_tigrfam = abundTable.transpose().dot(copyNumberTable_tigrfam_consolidated).transpose().round()
final_df_tigrfam = final_df_tigrfam.loc[final_df_tigrfam.sum(axis=1) != 0]
os.mkdir(cwd + "/TIGRFAM_metagenome")
final_df_tigrfam.to_csv(cwd + "/TIGRFAM_metagenome/TIGRFAM_metagenome.txt", sep="\t")


# annotate metagenomes with additional annotations
# KO
if (in_verbose):
    print("Anotating all predicted metagenomes")
final_df_annotated_ko = fp.addAnnotations(final_df, otherPath + "ko00001.txt")
final_df_annotated_ko.to_csv(cwd + "/KO_metagenome/KO_metagenome_with_description.txt", sep="\t")
# COG
final_df_annotated_cog = fp.addAnnotations(final_df_cog, otherPath + "cognames2003-2014.tab")
final_df_annotated_cog.to_csv(cwd + "/COG_metagenome/COG_metagenome_with_description.txt", sep="\t")
# Pfam
final_df_annotated_pfam = fp.addAnnotations(final_df_pfam, otherPath + "pfam_mapping.txt")
final_df_annotated_pfam.to_csv(cwd + "/Pfam_metagenome/Pfam_metagenome_with_description.txt", sep="\t")
# TIGRFAM
final_df_annotated_tigrfam = fp.addAnnotations(final_df_tigrfam, otherPath + "TIGRFAMs_9.0_INFO.txt")
final_df_annotated_tigrfam.to_csv(cwd + "/TIGRFAM_metagenome/TIGRFAM_metagenome_with_description.txt", sep="\t")

# predict kegg rections
df_kegg_reaction = fp.ec2Reaction_kegg(final_df_ec, otherPath + "kegg_ec2reaction")
os.mkdir(cwd + "/KO_metagenome/metabolome")
df_kegg_reaction.to_csv(cwd + "/KO_metagenome/metabolome/reaction_metagenome.txt", sep='\t')

# predict kegg metabolites
if(in_verbose):
    print("Predicting KEGG metabolic potential")
comp_refData = pd.read_csv(otherPath + "kegg_reaction2comp.gz", sep='\t', index_col=0, compression='gzip')
comp_refData = comp_refData.drop(comp_refData.columns.difference(final_df.index),axis=1)
final_df = final_df.drop(final_df.index.difference(comp_refData.columns))
df_comMetagenome = comp_refData.dot(final_df)

# add annotations
if(in_verbose):
    print("Annotating metabolites")
df_comMetagenome = fp.addAnnotations(df_comMetagenome, otherPath + 'kegg_compAnnotations.txt')
df_comMetagenome.to_csv(cwd + "/KO_metagenome/metabolome/comp_metagenome.txt", sep='\t')

# run MinPath
if (in_verbose):
    print("Running MinPath (KO)")
final_df_minPath = fp.runMinPath(final_df_annotated_ko, funpredPath, cwd + "/KO_metagenome", "kegg")
# groupby levels
fp.summarizeByFun(final_df_minPath, "A").to_csv(cwd + "/KO_metagenome/summarized_by_A.txt", sep="\t")
fp.summarizeByFun(final_df_minPath, "B").to_csv(cwd + "/KO_metagenome/summarized_by_B.txt", sep="\t")
fp.summarizeByFun(final_df_minPath, "C").to_csv(cwd + "/KO_metagenome/summarized_by_C.txt", sep="\t")
fp.summarizeByFun(final_df_minPath, "Pathway_Module").to_csv(cwd + "/KO_metagenome/summarized_by_Pathway_Module.txt",
                                                             sep="\t")

###################### MetaCyc ################3
if (in_verbose):
    print("Running MinPath (RXN)")
final_df_ec_rxn = fp.ec2RXN(final_df_ec, otherPath + "ec2rxn_new")
final_df_ec_rxn.to_csv(cwd + "/MetaCyc_metagenome/RXN_metagenome.txt", sep="\t")
os.mkdir(cwd + "/MetaCyc_metagenome" + "/minPath_files")

#run MinPath EC
final_PathwayAbundance_df = fp.runMinPath(final_df_ec_rxn, funpredPath, cwd + "/MetaCyc_metagenome/minPath_files",
                                          "metacyc")
final_PathwayAbundance_df.to_csv(cwd + "/MetaCyc_metagenome/PathwayAbundance.table", sep="\t")
final_PathwayAbundance_df = fp.addMetaCycPathwayName(final_PathwayAbundance_df, otherPath + "path_to_Name.txt")
final_PathwayAbundance_df.to_csv(cwd + "/MetaCyc_metagenome/PathwayAbundance_with_names.table", sep="\t")
fp.summarizeByFun(final_PathwayAbundance_df, "Type").to_csv(
    cwd + "/MetaCyc_metagenome/Pathway_summarize_by_Types.table", sep="\t")

######BacDive Phenotype
abundTable = pd.read_csv(in_otuTab,sep="\t",index_col=0)
if(in_phenotype):
    if(in_verbose):
        print("Predicting phenotypic data of OTUs/ASVs.")
    os.mkdir(cwd + "/PhenotypeData")
    fp.blast(in_repSet, dbPath + "blastDB/16S_LPSN", in_numCores, 100, cwd + '/PhenotypeData')
    df_pheno_blastOut = pd.read_csv(cwd+'/PhenotypeData/blast.outfmt6',sep='\t',header=None,index_col=0)
    #subset bacdive
    df_bacdive_sub = fp.subsetBacdive(df_pheno_blastOut,abundTable,dataPath +'phenotype/bacdive_newest_edited.txt.gz')
    df_bacdive_sub.transpose().to_csv(cwd+'/PhenotypeData/bacDive_phenoData.txt',sep='\t')
timeTaken = round(time.time() - start_time, 2)

print("Completed pipline in " + str(timeTaken) + " seconds.")
"""
