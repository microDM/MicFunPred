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
# warnings.filterwarnings("ignore")
import os
import shutil
import time
import subprocess
import numpy as np
import random

# make options
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--otu_table", metavar="PATH",type=str,
                  help="[Required] Tab delimited OTU table")
parser.add_argument("-r", "--repset_seq",metavar="PATH",type=str,
                  help="[Required] Multi-fasta file of OTU/ASV sequences")
parser.add_argument("-p", "--perc_ident",type=str,
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

################################ Definition ################################3

def makeKOTabke(df,abundData,coreNum,percentGenomes):
    """
    Consolidate the KO copy number table with respect to OTU table
    Args:
        df: Gene copy number table (dataframe)
        abundData: OTU abundance table (output of selectBlastHits_assignGenus_subsetOtuTable)
        coreNum: value in range 0 to 1. If a gene present in coreNum*genus, then it will be considered as core gene.

    Returns:
        DataFrame: Consolidated gene copy number table (dataframe)
    """
    taxonomyList = list(abundData.index)
    dfToReturn = pd.DataFrame()
    for tax in taxonomyList:
        temp_df = df[df.index.str.contains(tax,na=False)]
        if(int(round(temp_df.shape[0]*percentGenomes)) < 1):
            k = 1
        else:
            k = int(round(temp_df.shape[0]*percentGenomes))
        temp_df = temp_df.loc[random.sample(list(temp_df.index),k=k)]
        numberOfKO = 0
        keepCols = list()
        colsToDrop = list()
        for col in temp_df.columns:
            if(coreNum == 0):
                if(np.count_nonzero(temp_df[col]) >= 1):
                    keepCols.append(col)
                    numberOfKO = numberOfKO + 1
                else:
                    colsToDrop.append(col)    
            else:
                if(np.count_nonzero(temp_df[col]) >= round(temp_df.shape[0]*coreNum)):
                    keepCols.append(col)
                    numberOfKO = numberOfKO + 1
                else:
                    colsToDrop.append(col)
        #retain columns from keepCols
        temp_df = temp_df.drop(colsToDrop,axis=1)
        median_series = temp_df.mean()
        median_series[median_series.between(0,1,False)] = 1
        #median_df = pd.Series.to_frame(median_series).transpose().round()
        median_df = pd.Series.to_frame(median_series).transpose()
        dfToReturn = dfToReturn.append(median_df, ignore_index = True,sort=False)
    dfToReturn.index = taxonomyList
    #replace NA with 0
    dfToReturn = dfToReturn.fillna(0)
    return dfToReturn

###########################################        MAIN         ######################33

#data path
funpredPath = fp.baseDir()

cwd = os.getcwd()
dbPath = funpredPath + "/data/db/"
dataPath = funpredPath + "/data/"
otherPath = funpredPath + "/data/other/"
#create databases
if(not os.path.exists(funpredPath+'/data/db/blastDB')):
    print("Please wait. Creating databases for first time. This is one-time process.")
    os.chdir(funpredPath + '/data/db')
    os.mkdir('blastDB')
    files = [x for x in os.listdir(funpredPath+'/data/db') if '.fasta.gz' in x]
    for i in files:
        if('all' in i):
            dbname = '16S_all'
        else:
            dbname = '16S_LPSN'
        cmd = 'gunzip -c ' + i + ' | makeblastdb -dbtype nucl -input_type fasta -in - -title 16Sdb -out blastDB/' + dbname
        print(cmd)
        subMakeBlastdb = subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE,shell=True)
        subMakeBlastdb.wait()
else:
    pass

os.chdir(cwd)

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
    fp.blast(in_repSet, dbPath + "blastDB/16S_all", in_numCores, in_percIdentCutOff, cwd)
    tax_dict, abundTable = fp.selectBlastHits_assignGenus_subsetOtuTable(cwd + "/blast.outfmt6", in_otuTab,
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

#percent prediction
os.mkdir(cwd + "/KO_metagenome")
for p in range(10,110,10):
    print(p)
    copyNumberTable_ko_consolidated = makeKOTabke(copyNumberTable_ko, abundTable, in_coreperc,p*0.01).fillna(0)
    copyNumberTable_ko_consolidated.to_csv(cwd + "/predicted_KO_" + str(p) + ".txt", sep="\t")
    if (in_verbose):
        print("Predicting KO, EC, Pfam, COG and TIGRFAM metagenome")
    final_df = abundTable.transpose().dot(copyNumberTable_ko_consolidated).transpose()
    final_df = final_df.loc[final_df.sum(axis=1) != 0]
    final_df.to_csv(cwd + "/KO_metagenome/KO_metagenome_" + str(p) + ".txt", sep="\t")