#!/usr/bin/python
import plotly.graph_objects as go
import pandas as pd
from optparse import OptionParser
import colorlover as cl
import plotly_express as px

# add options
parser = OptionParser()
parser.add_option('-i', '--input-file', help='[Required] Tab delimited abundace file (Ex: KO_metagenome.txt',
                  dest='in_inFile')
parser.add_option('-n', '--num-samples',help='[Required] Number of samples in input-file',
                  dest='in_numSamples')
parser.add_option('-l', '--level', help='[Required] Annotation level to be plotted',
                  dest='in_level')
parser.add_option('-m', '--map-file', help='[Required] Tab delimited metadata file',
                  dest='in_mapFile',default=None)
parser.add_option('-v', '--map-variable', help='[Required] metadata varibale name. If not given, then 1st '
                  'column from map-file (excluding sample name) will be considered.',
                  dest='in_metadatavar', default=None)

(options,args) = parser.parse_args()
requiredOptions = 'in_inFile in_numSamples in_level in_mapFile in_metadatavar'.split()
for r in requiredOptions:
    if(options.__dict__[r] is None):
        parser.error('Please give all required arguments. %s is missing'%r.upper())
        exit()
in_inFile = options.in_inFile
in_numSamples = int(options.in_numSamples)
in_level = options.in_level
in_mapFile = options.in_mapFile
in_metadataVar = options.in_metadatavar


def prepareDataframeDict(mainDf,mapDf,metadataVar,level,numOfSamples):
    '''
    :param mainDf: Metagenome dataframe with levels
    :param mapDf: mapping dataframce
    :param metadataVar: metadata variable to melt
    :param level: level name in mainDf
    :param numOfSamples: number of samples
    :return: dictionary of dataframes
    '''
    df_dict = dict()
    temp = mainDf[mainDf.columns[0:numOfSamples]]
    temp['level'] = mainDf[level]
    temp = temp.groupby('level').sum()
    temp = temp/temp.sum()
    temp = temp.transpose().join(mapDf).drop(mapDf.columns.difference([metadataVar]),axis=1).melt(metadataVar)
    temp.columns = [metadataVar,level,'Relative_abundance']
    buttonNames = list(set(temp[level]))
    #print(buttonNames)
    for i in set(temp[level]):
        #print(i)
        #print(temp[temp[level].str.contains(i)])
        df_dict[i] = temp[temp[level] == i]
        print(i)
    return df_dict,buttonNames

def addTrace(df_dict,metadataVar,level):
    '''
    :param df_dict: Dictionary of dataframes
    :return: plotly.graphics_obj with added traces
    '''
    c = cl.scales['12']['qual']['Paired']
    fig = go.Figure()
    for k,v in sorted(df_dict.items()):
        df_temp = df_dict[k]
        x = df_temp[metadataVar].to_list()
        y = df_temp['Relative_abundance'].to_list()
        name = k
        fig.add_trace(go.Box(x=x,y=y,visible=False,))
    #print(fig.data)
    return fig

def prepareAnnotations(df_dict,metadatavar,level):
    annotations = dict()
    for k, v in sorted(df_dict.items()):
        df_temp = df_dict[k]
        temp = []
        temp.append(
            dict(
                x=df_temp[metadatavar].to_list(),
                y=df_temp['Relative_abundance'].to_list(),
                text="Relative Abundance of %s level"%k
            )
        )
        annotations[k] = temp
    return annotations

def prepareButtons(annotations,buttonNames):
    buttons = []
    for j,i in enumerate(sorted(buttonNames)):
        # make booleanList
        tempBooleanList = [False for i in buttonNames]
        tempBooleanList[j] = True
        temp = dict(
            label=i,
            method='update',
            args=[{'visible':tempBooleanList},
                  {'title':'Relative abundance of %s'%i},
                  {'annotations':annotations[i]}]
        )
        buttons.append(temp)
    return buttons

def addUpdateMenus(fig,annotations,buttons):
    fig.update_layout(
        updatemenus=[
            go.layout.Updatemenu(
                buttons=buttons,
                direction="down",
                xanchor='left',
                yanchor='bottom'
            )]
    )
    return fig

df = pd.read_csv(in_inFile,sep='\t',index_col=0)
mapdf = pd.read_csv(in_mapFile,sep='\t',index_col=0)
df_dict,buttonNames = prepareDataframeDict(df,mapdf,in_metadataVar,in_level,12)
fig = addTrace(df_dict,in_metadataVar,in_level)
annotations = prepareAnnotations(df_dict,in_metadataVar,in_level)

buttons = prepareButtons(annotations,buttonNames)

fig.update_layout(
        updatemenus=[
            go.layout.Updatemenu(
                buttons=buttons,
            )]
    )
fig.show()
