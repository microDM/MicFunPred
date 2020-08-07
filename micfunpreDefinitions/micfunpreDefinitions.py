
'''
Created on 22-Feb-2019

@author: Dattatray Mongad,
                    Junior Research Fellow (PhD Student),
                    National Centre fo Microbial Resource,
                    National Centre for Cell Science,
                    Pune, Maharshtra, India.
'''

import pandas as pd
import numpy as np
import subprocess
import os
import re

def baseDir():
    """Return baseDir of installation
    """
    return os.path.dirname(os.path.abspath(__file__))

def blast(query,db,in_numCores,cwd):
    """
    Do BLAST with provided query and database
    Args:
        query: Query file
        db: db path/name for doing blast
        in_numCores: Number of cores to use for BLAST
        in_percIdentCutOff: Percent identity cutoff
        cwd: Current working directory

    Returns:
        None
    """
    cmd = "blastn -out " + cwd + "/out.blast -outfmt 6 -query " + query + " -db " + db + " -num_threads " + in_numCores + " -max_target_seqs 1"
    subBlast = subprocess.Popen(cmd,stdin=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    output,error = subBlast.communicate()
    return None

def selectBlastHits_assignGenus_subsetOtuTable(blastOut,otuTable,blastcutoff,level):
    """
    Select best hits from BLAST output and assign taxonmy with respect to given identity cutoff
    Args:
        blastOut: Filtered blast output file (outfmt 6)
        otuTable: OTU table
        blastcutoff: percent identity cutoff to assign genus to the sequences

    Returns:
        DataFrame: OTU/ASV table containg OTUs/ASVs which have assigned taxonomy through BLAST
    """
    df = pd.read_csv(blastOut,sep='\t',index_col=0,header=None)
    df = df.loc[df[2]>=blastcutoff]
    if(df.empty):
        exit('None of OTU/ASV sequences passed the given percent identity cut-off')
    otuId_tax_dict = dict(zip(df.index,df[1].str.split("_",expand=True)[1]))
    df = pd.read_csv(otuTable,sep="\t",index_col=0)
    df.columns = df.columns.astype(str)
    df = df.reindex(list(otuId_tax_dict.keys()))
    df["taxonomy"] = list(otuId_tax_dict.values())
    df = df.groupby(["taxonomy"]).sum()
    return otuId_tax_dict,df
     

def makeTable16S(df,func,taxonomyList):
    """
    Consolidate the 16S rRNA copy number table with respect to taxonomy
    Args:
        df: Dataframe of 16S rRNA copy numner
        func: mean/mode/median to be taken
        taxonomyList: list of taxonomy

    Returns:
        DataFrame: consolidated dataframe for organisms provided in taxonomy list
    """
    #make table
    temp_dict = dict()
    for tax in taxonomyList:
        num = round(float(df[df.index.str.contains(tax,na=False)].mean()),2)
        if(num == 0):
            temp_dict[tax] = 1
        else:
            temp_dict[tax] = num
    df_consolidated = pd.DataFrame.from_dict(temp_dict,orient="index")
    del(temp_dict)
    return df_consolidated

def makeKOTable(df,abundData,coreNum):
    """
    Consolidate the KO copy number table with respect to OTU table
    Args:
        df: Gene copy number table (dataframe)
        abundData: OTU abundance table (output of selectBlastHits_assignGenus_subsetOtuTable)
        coreNum: value in range 0 to 1. If a gene present in coreNum*genus, then it will be considered as core gene.

    Returns:
        DataFrame: Consolidated gene copy number table (dataframe)
    """
    # old
    taxonomyList = list(abundData.index)
    dfToReturn = pd.DataFrame()
    for tax in taxonomyList:
        temp_df = df[df.index.str.contains(tax,na=False)]
        n = round(temp_df.shape[0]*coreNum)
        temp_df = temp_df[temp_df.columns[temp_df.astype(bool).sum()>=n]]
        median_series = temp_df.mean()
        median_series[median_series.between(0,1,False)] = 1
        #median_df = pd.Series.to_frame(median_series).transpose().round()
        median_df = pd.Series.to_frame(median_series).transpose()
        dfToReturn = dfToReturn.append(median_df, ignore_index = True,sort=False)
    dfToReturn.index = taxonomyList
    #replace NA with 0
    dfToReturn = dfToReturn.fillna(0)
    """
    # new
    # predict cut-offs
    taxonomyList = list(abundData.index)
    dfToReturn = pd.DataFrame()
    cutOffDict = dict()
    for tax in taxonomyList:
        half = df_cutOff.loc[tax]*coreNum
        cutOff = int(df_cutOff.loc[tax].lt(half).idxmax())*0.01
        cutOffDict[tax] = cutOff
    # predict genes based on predicted cut-offs
    for tax in taxonomyList:
        temp_df = df[df.index.str.contains(tax,na=False)]
        n = round(temp_df.shape[0]*cutOffDict[tax])
        temp_df = temp_df[temp_df.columns[temp_df.astype(bool).sum()>=n]]
        median_series = temp_df.median()
        median_series[median_series.between(0,1,False)] = 1
        #median_df = pd.Series.to_frame(median_series).transpose().round()
        median_df = pd.Series.to_frame(median_series).transpose()
        dfToReturn = dfToReturn.append(median_df, ignore_index = True,sort=False)
    dfToReturn.index = taxonomyList
    #replace NA with 0
    dfToReturn = dfToReturn.fillna(0)
    """
    return dfToReturn

def addAnnotations(metagenomeDf,keggFile):
    """
    Add Kegg Annotations to predicted metagenome profile
    Args:
        metagenomeDf: predicted metagenome profile DataFrame

    Returns: Predicted metagenome with KEGG annotations

    """
    # read kegg annotations
    kodf = pd.read_csv(keggFile, sep="\t", index_col=0,engine='python')
    metagenomeDf = metagenomeDf.join(kodf)
    return metagenomeDf

def summarizeByFun(metagenomeDf,group):
    """
    Consolidate on the basis of group
    Args:
        metagenomeDf: Annotated Metagenme matrix
        group: string by which dataFrame is to be categorized

    Returns: Consolidated DataFrame

    """
    return metagenomeDf.groupby(group).sum()

def runMinPath(metagenomeDf,funpredPath,outPath,typeOfPrediction):
    """
    Run MinPath
    Args:
        metagenomeDf: Predicted Meatagenome matrix
        funpredPath: Path of funPred
        outPath: path to store files
        typeOfPrediction: kegg or metacyc depending on type of input KO or EC

    Returns:
        DataFrame: Pruned metagenome content dataframe
    """
    #make input for minPath and run MinPath
    if(typeOfPrediction == "kegg"):
        with open(outPath + "/minpath_in.ko","w") as f:
            for i,j in enumerate(list(metagenomeDf.index)):
                    f.write(str(i) + "\t" + str(j)+"\n")
        cmd = "python2 " + funpredPath  +  "/MinPath1.4_micfunpred.py " + funpredPath + " " + outPath + " -ko " + outPath + "/minpath_in.ko -report " + outPath + "/minpath.out"
        a = os.popen(cmd).read()
        minPathOutDF = pd.read_csv(outPath + "/minpath.out", sep="\t", index_col=0, header=None)
        predictedMaps = ["ko" + x.replace("path ", "") for x in
                         list(minPathOutDF[minPathOutDF[3] == "minpath 1"].index)]
        temp = pd.DataFrame()
        for i in predictedMaps:
            temp = temp.append(metagenomeDf[metagenomeDf["C"].str.contains(i, na=False)])
        temp.to_csv(outPath + "/KO_metagenome_minPath_pruned.txt", sep="\t")
        #os.remove(outPath + "/temp.ko")
        return temp
    elif(typeOfPrediction == "metacyc"):
        minpathOutFiles = []
        for sampleName in metagenomeDf.columns:
            minPtahInFile = outPath + "/" + sampleName + "_minpath_in.txt"
            #make list of files
            minpathOutFiles.append(outPath + "/" +sampleName + "_minpath.out.details")
            #create input file and run MinPath
            minPtahInFile_fh = open(minPtahInFile,"w")
            minPtahInFile_fh.writelines(['read' + str(i) + "\t" + str(j) + "\n" for i,j in enumerate(list(metagenomeDf[metagenomeDf[sampleName]>0].index))])
            cmd = "python2 " + funpredPath + "/MinPath1.4_micfunpred.py " + funpredPath + " " + outPath + " -any " + minPtahInFile + " -map " + funpredPath + "/data/path_to_RXN.txt -report " + outPath + "/" + sampleName + "_minpath.out -details " + outPath + "/" + sampleName + "_minpath.out.details"
            a = os.popen(cmd).read()
        #create pathway abundance dataframe from all files
        metagenomeDf_reindexed = metagenomeDf.copy()
        pathDict = {}
        for i in minpathOutFiles:
            path = ""
            rxn = ""
            iFH = open(i,"r")
            for line in iFH.readlines():
                matchObj = re.match("^path.*\#\s(\S+)",line)
                if(matchObj):
                    path = matchObj.group(1)
                matchObj = re.match("^\s+(\S+)",line)
                if(matchObj):
                    rxn = matchObj.group(1)
                    pathDict[rxn] = path
            iFH.close()
        # select rows having RXN found
        metagenomeDf_reindexed = metagenomeDf_reindexed.loc[list(pathDict.keys())]
        metagenomeDf_reindexed.index = metagenomeDf_reindexed.index.to_series().replace(pathDict)
        return metagenomeDf_reindexed.groupby(metagenomeDf_reindexed.index).sum()
    
def ec2RXN(df,ec2rxnFile):
    """

    Args:
        df: EC metagenome dataframe
        ec2rxnFile: map file for ec to RXN

    Returns:

    """
    # make ec dict
    ec2RXNDict = dict()
    with open(ec2rxnFile,"r") as f:
        for i in f.readlines():
            if("-" in i.split("\t")[0]):
                ec2RXNDict[i.split("\t")[0][:-2]] = i.split("\t")[1].strip()
            else:
                ec2RXNDict[i.split("\t")[0]] = i.split("\t")[1].strip()
    #iterate over metagneome dataframe
    tempDf = pd.DataFrame(columns=list(df.columns))
    for ec in df.index:
        if("-" in ec):
            ecTemp = ec[:-2]
        else:
            ecTemp = ec
        if(ecTemp in ec2RXNDict.keys()):
            rxnList = ec2RXNDict[ecTemp].split(",")
            for rxn in ec2RXNDict[ecTemp].split(","):
                tempDf.loc[rxn] = list(df.loc[ec])
    return tempDf

def addMetaCycPathwayName(pathAbundanceDf,pathNameFile):
    """
    :param pathAbundanceDf: Pathway abundance dataframe
    :param pathNameFile: Reference File of pathways
    :return: Pathway abundance dataframes with pathway types and common-names
    """
    pathNameDf = pd.read_csv(pathNameFile,sep="\t",index_col=0)
    return pathAbundanceDf.join(pathNameDf)

def ec2Reaction_kegg(df,ec2reactionFile):
    # make ec dict
    ec2ReactionDict = dict()
    with open(ec2reactionFile,"r") as f:
        for i in f.readlines():
            if("-" in i.split("\t")[0]):
                ec2ReactionDict[i.split("\t")[0][:-2]] = i.split("\t")[1].strip()
            else:
                ec2ReactionDict[i.split("\t")[0]] = i.split("\t")[1].strip()
    #iterate over metagneome dataframe
    tempDf = pd.DataFrame(columns=list(df.columns))
    for ec in df.index:
        if("-" in ec):
            ecTemp = ec[:-2]
        else:
            ecTemp = ec
        if(ecTemp in ec2ReactionDict.keys()):
            rxnList = ec2ReactionDict[ecTemp].split(",")
            for rxn in ec2ReactionDict[ecTemp].split(","):
                tempDf.loc[rxn] = list(df.loc[ec])
    return tempDf

def subsetBacdive(df_blast,abundanceTable,bacdiveFileName):
    """
    :param df_blast:
    :param bacdiveFileName:
    :param abundanceTable
    :return:
    """
    #relative abundance in abundance table
    abundanceTable = abundanceTable.div(abundanceTable.sum())
    abundanceTable['ASV'] = abundanceTable.index
    #read bacdive file
    df_bacdive = pd.read_csv(bacdiveFileName,sep='\t',index_col=0,compression='gzip',low_memory=False).transpose()
    found_acc_num = list(df_blast[1].str.replace('\.\d+',''))
    dfList = []
    tempList = []
    for i,j in zip(found_acc_num,df_blast.index):
        if( not df_bacdive[df_bacdive['all_reference_acc_num'].str.contains(i, na=False)].empty):
            tempDf = df_bacdive[df_bacdive['all_reference_acc_num'].str.contains(i, na=False)]
            dfList.append(tempDf)
            for i in range(tempDf.shape[0]):
                tempList.append(j)
    df_bacdive = pd.concat(dfList,sort=True)
    df_bacdive['bacdiveID'] = df_bacdive.index
    df_bacdive['seqName'] = tempList
    df_bacdive.dropna(axis=1,how='all',inplace=True)
    df_bacdive.index.name = 'bacdiveID'
    #add abundance values
    df_bacdive = df_bacdive.merge(abundanceTable, left_on='seqName', right_on='ASV').drop('seqName', axis=1).set_index('ASV')
    return df_bacdive