#!/usr/bin/python

"""
Created on 08-Aug-2019

@author: Dattatray Mongad,
                    Junior Research Fellow (PhD Student),
                    National Centre fo Microbial Resource,
                    National Centre for Cell Science,
                    Pune, Maharshtra, India.
"""

import pandas as pd
import plotly as py
import plotly_express as px
from optparse import OptionParser
import os

# Add options

parser = OptionParser()
parser.add_option('-i','--input',dest='in_inputData',
                  help='[Required] Tab deliminated input data. Example: Output created by MicFunpred '
                       '(KO_metagenome/KO_metagenome_with_description.txt')
parser.add_option('-l','--level',dest='in_levels',help='[Required] Level in input data to collapse with. '
                  'Multiple levels can be given in comma-separated format. Example: '
                  'A,B,C')
parser.add_option('-m','--map',dest='in_mapData',default='',
                  help='Tab deliminated metadata. First column should represent '
                  'same column names present in input data. If not provided plots x-axis will have sample names.')
parser.add_option('-v','--metadata-variable',default='',dest='in_metadataVariable',help='Variable name '
                  'from metadata/map file to be used for grouping samples. Multiple variables '
                  'can be given in comma-separated format. Do not give'
                    ' if metadata not provided. Example: group1,group2,group3')
parser.add_option('-n','--num-samples',dest='in_numSamples',default=1,help='Number of samples '
                  'present in input data. If map data is not passed [Required], it will be used. If map data '
                  'is passed, not necessary to give this argument.')
parser.add_option('-o','--out-dir',dest='in_outDir',help='Output directory.',default='plotBar_out')

# prepare options
(options, args ) = parser.parse_args()
required = 'in_inputData in_mapData in_metadataVariable in_levels'.split()
for i in required:
    if options.__dict__[i] is None:
        parser.error('Argument %s is required'%i)
    else:
        in_mapData = options.in_mapData
        in_metadataVariable = options.in_metadataVariable.split(',')
        in_inputData = options.in_inputData
        in_levels = options.in_levels.split(',')
        in_outDir = options.in_outDir
        in_numSamples = int(options.in_numSamples)

def prepareDataFrame(mainDf,mapDf,levelName,metadataVariable,numSamples):
    '''

    :param df: Annotated dataframe
    :param metadataVariable: Varibale name to plot
    :return: melted dataframe having levelName and metadataVariable
    '''
    if(mapDf == '' and metadataVariable == ''):
        mainDf['temp'] = mainDf.index
        mainDf = mainDf.groupby(['temp',levelName]).sum().reset_index().groupby(levelName).sum()
        mainDf = mainDf.div(mainDf.sum())
        mainDf = mainDf[list(mainDf.columns)[:in_numSamples]]
        mainDf[levelName] = mainDf.index
        mainDf = mainDf.melt(levelName)
        mainDf.columns = [levelName,'SampleName','value']
    else:
        mainDf = mainDf.groupby(levelName).sum()  # groupby level
        # add maetadata and groupby according to it
        mainDf['temp'] = mainDf.index
        mainDf = mainDf.groupby(['temp',levelName]).sum().reset_index().groupby(levelName).sum()
        mainDf = mainDf.div(mainDf.sum())
        mainDf = mainDf.transpose().join(mapDf)
        mainDf = mainDf.drop(mapDf.columns.difference([metadataVariable]),axis=1)
        mainDf = mainDf.groupby(metadataVariable).sum()
        mainDf[metadataVariable] = mainDf.index
        mainDf = mainDf.melt(metadataVariable)
        mainDf.columns = [metadataVariable,levelName,'Relative_proportion']
    return mainDf

# make output direcory
if(os.path.exists(in_outDir)):
    os.system('rm -r ' + in_outDir)
    os.mkdir(in_outDir)
else:
    os.mkdir(in_outDir)

# read data
df = pd.read_csv(in_inputData,sep='\t',index_col=0)
if(in_mapData != ''):
    df_map = pd.read_csv(in_mapData,sep='\t',index_col=0)

if(in_mapData or in_metadataVariable is None):
    for m in in_metadataVariable:
        for l in in_levels:
            temp = in_inputData.split('.')[0]
            outFileName = in_outDir + '/' + temp + '_' + m + '_' + l + '.html'
            dftemp = prepareDataFrame(df,df_map,l,m,'')
            fig = px.bar(data_frame=dftemp,x=m,y='Relative_proportion',color=l)
            with open(outFileName,'w') as f:
                f.write(fig.to_html())
elif(not ((in_mapData or in_metadataVariable) is None) ):
    for l in in_levels:
        temp = in_inputData.split('.')[0].split('/')[-1]
        outFileName = in_outDir + '/' + temp + '_' + l + '.html'
        print(temp)
        dftemp = prepareDataFrame(df,'',l,'',in_numSamples)
        fig = px.bar(data_frame=dftemp,x='SampleName',y='value',color=l,width=1000)
        with open(outFileName,'w') as f:
            f.write(fig.to_html())
