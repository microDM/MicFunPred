#!/usr/bin/python
import plotly.graph_objects as go
import pandas as pd
from optparse import OptionParser
import os

# add options
parser = OptionParser()
parser.add_option('-i', '--input-file',
                  help='[Required] Tab delimited abundace file (Ex: KO_metagenome.txt', dest='in_inFile')
parser.add_option('-n', '--num-samples', help='[Required] Number of samples in input-file',
                  dest='in_numSamples')
parser.add_option('-m', '--map-file', help='[Optional] Tab delimited metadata file',
                  dest='in_mapFile', default=None)
parser.add_option('-v', '--map-variable', help='[Optional] metadata varibale name. If not given, then 1st '
                  'column from map-file (excluding sample name) will be considered.',
                  dest='in_metadatavar')
parser.add_option('-l', '--levels', help='[Optional] Comma separated annotation levels to be plotted',
                  dest='in_levels')
(options, args) = parser.parse_args()
# check required options
requiredOptions = 'in_inFile in_numSamples'.split()
for r in requiredOptions:
    if(options.__dict__[r] is None):
        parser.error('Please give all required arguments. %s is missing' % r.upper())
        exit()
in_inFile = options.in_inFile
in_numSamples = int(options.in_numSamples)
in_mapFile = options.in_mapFile
if(not options.__dict__['in_levels'] is None):
    in_levels = options.in_levels.split(',')
else:
    in_levels = options.in_levels
if(not options.__dict__['in_metadatavar'] is None):
    in_metadatavar = options.in_metadatavar
else:
    in_metadatavar = None

# definitions #####################3


def prepareDataframeDict(maindf, mapdf, metadavar, levels, num):
    '''
    :param maindf: Dataframe of metagenome
    :param mapdf: Mapping dataframe
    :param metadavar: Variable from mapping dataframe
    :param levelName: Annotation level name or list of levels
    :param num: number of samples
    :return: dictionary of melted dataframes with keys as level names and values as melted dataframes
    '''
    df_dict = {}
    if(mapdf.empty and metadavar == None):
        for l in levels:
            temp = maindf.copy()
            temp['index'] = temp.index
            temp = temp.groupby([l, 'index'], axis=0).sum().reset_index().groupby(l).sum()
            temp = temp/temp.sum()
            temp['index'] = temp.index
            temp = temp.melt('index')
            temp.columns = [l, 'Samples', 'Relative_Abundance']
            df_dict[l] = temp
        return df_dict
    else:
        for l in levels:
            temp = maindf.copy()
            temp['index'] = temp.index
            temp = temp.groupby([l, 'index'], axis=0).sum().reset_index().groupby(l).sum()
            temp = temp/temp.sum()
            temp = temp.transpose().join(mapdf)
            temp = temp.drop(mapdf.columns.difference([metadavar]), axis=1).melt(
                metadavar).groupby(['variable', metadavar]).sum().reset_index()
            temp.columns = [l, 'Samples', 'Relative_Abundance']
            # print(temp['Samples'])
            df_dict[l] = temp
        return df_dict


def addTrace(df_dict):
    '''
    :param df_dict: Dictionary of dataframes to be plotted
    :return: Plotly.graph_object with traces
    '''
    fig = go.Figure()
    for k, v in sorted(df_dict.items()):
        df_temp = df_dict[k]
        nameList = list(set(df_temp[k]))
        for i, j in enumerate(nameList):
            x = df_temp[df_temp[k] == j]['Samples'].to_list()
            y = df_temp[df_temp[k] == j]['Relative_Abundance'].to_list()
            name = nameList[i]
            fig.add_trace(go.Bar(x=x, y=y, name=name, visible=False))
    fig.update_layout(barmode='stack')
    return fig


def prepareAnnotations(df_dict):
    '''
    :param df_dict:  Dictionary of dataframes to be plotted
    :return:  Dictionary of annotations
    '''
    annotations = dict()
    for k, v in sorted(df_dict.items()):
        df_temp = df_dict[k]
        temp = []
        for i, j in enumerate(set(df_temp[k])):
            temp.append(
                dict(
                    x=df_temp[df_temp[k] == j]['Samples'].to_list(),
                    y=df_temp[df_temp[k] == j]['Relative_Abundance'].to_list(),
                    text="Relative Abundance"
                )
            )
        annotations[k] = temp
    return annotations


def prepareButtons(annotations, buttonNames):
    '''
    :param annotations: Dictionary of annotaions
    :return: List of buttons
    '''
    buttons = []
    # making boolean list
    lengthDict = {k: len(value) for k, value in annotations.items()}
    hugeBooleanList = []
    for k, v in lengthDict.items():
        hugeBooleanList = hugeBooleanList + [False for i in range(lengthDict[k])]
    start = 0
    end = 0
    for i, j in enumerate(buttonNames):
        tempBooleanList = hugeBooleanList.copy()
        end = start + lengthDict[j]
        tempBooleanList[start:end] = [True for i in range(end - start)]
        temp = dict(
            label=j,
            method='update',
            args=[{'visible': tempBooleanList},
                  {'title': 'Relative Abundance of %s level' % j,
                   'annotations': annotations[j]}]
        )
        buttons.append(temp)
        start = end
    return buttons


def addUpdateMenus(fig, annotations, buttons):
    fig.update_layout(
        updatemenus=[
            go.layout.Updatemenu(
                buttons=buttons,
                direction="down",
                xanchor='left',
                yanchor='bottom'
            )],
        margin=go.layout.Margin(
            l=50,
            r=10,
            pad=4
        )
    )
    return fig


###################### MAIN #########################
df = pd.read_csv(in_inFile, sep='\t', index_col=0)  # read input dataframe
# read levels
if(in_levels is None):
    levels = sorted(list(df.columns)[in_numSamples:])

# prepare dataframe (data to plot)
if(in_mapFile is None):
    df_dict = prepareDataframeDict(df, pd.DataFrame(), None, levels, in_numSamples)
    fig = addTrace(df_dict)
    annotations = prepareAnnotations(df_dict)
    buttons = prepareButtons(annotations, levels)
    fig = addUpdateMenus(fig, annotations, buttons)
    fig.show()

if(in_mapFile is not None):
    mapdf = pd.read_csv(in_mapFile, sep='\t', index_col=0)
    if(in_metadatavar is None):
        in_metadatavar = list(mapdf.columns)[0]
    if(len(in_metadatavar.split(',')) == 1):
        df_dict = prepareDataframeDict(df, mapdf, in_metadatavar, levels, in_numSamples)
        fig = addTrace(df_dict)
        annotations = prepareAnnotations(df_dict)
        buttons = prepareButtons(annotations, levels)
        fig = addUpdateMenus(fig, annotations, buttons)
        fig.show()
    if(len(in_metadatavar.split(',')) > 1):
        if(os.path.exists('plots')):
            os.system('rm -r plots')
        os.mkdir('plots')
        for m in in_metadatavar.split(','):
            df_dict = prepareDataframeDict(df, mapdf, m, levels, in_numSamples)
            fig = addTrace(df_dict)
            annotations = prepareAnnotations(df_dict)
            buttons = prepareButtons(annotations, levels)
            fig = addUpdateMenus(fig, annotations, buttons)
            outputFile = 'plot_' + m + '.html'
            with open('plots/'+outputFile, 'w') as f:
                f.write(fig.to_html())
"""
print('Preparing dataframes')
#making dict of melted dataframes
if(in_mapFile is None):
    for i in metaDataVariables:
        temp = df.copy()
        temp['index'] = temp.index
        #consolidate
        temp = temp.groupby([i,'index'],axis=0).sum().reset_index().groupby(i).sum()
        #convert to relative abundance
        temp = temp/temp.sum()
        temp['index'] = temp.index
        temp = temp.melt('index')
        temp.columns = [i,'Samples','Relative_Abundance']
        df_dict[i] = temp
else:
    mdf = pd.read_csv(in_mapFile,sep='\t',index_col=0)

print('Adding trace')
fig = go.Figure()
#adding traces
for k,v in sorted(df_dict.items()):
    df_temp = df_dict[k]
    nameList = list(set(df_temp[k]))
    for i,j in enumerate(nameList):
        x = df_temp[df_temp[k] == j]['Samples'].to_list()
        y = df_temp[df_temp[k] == j]['Relative_Abundance'].to_list()
        name = nameList[i]
        fig.add_trace(go.Bar(x=x, y=y, name=name, visible=False))
fig.update_layout(barmode='stack')

print('Making annotations')
#making annoations
annotations = dict()
for k,v in sorted(df_dict.items()):
    df_temp = df_dict[k]
    temp = []
    for i,j in enumerate(set(df_temp[k])):
        temp.append(
            dict(
                x = df_temp[df_temp[k] == j]['Samples'].to_list(),
                y = df_temp[df_temp[k] == j]['Relative_Abundance'].to_list(),
                text = "Relative Abundance"
            )
        )
    annotations[k] = temp

#buttons
print('Preparing buttons')
buttons = []
#making boolean list
lengthDict = {k:len(value) for k,value in annotations.items()}
print(lengthDict)
hugeBooleanList = []
for k,v in lengthDict.items():
    hugeBooleanList = hugeBooleanList + [False for i in range(lengthDict[k])]

start = 0
end = 0
for i,j in enumerate(metaDataVariables):
    tempBooleanList = hugeBooleanList.copy()
    end = start + lengthDict[j]
    tempBooleanList[start:end] = [True for i in range(end-start)]
    print(j,start,end,str(tempBooleanList.count(True)))
    temp = dict(
        label=j,
        method='update',
        args=[{'visible': tempBooleanList},
              {'title': 'Abundance',
               'annotations': annotations[j]}]
    )
    buttons.append(temp)
    start = end + 1

#update layouts
print('Updating layout')
fig.update_layout(
    updatemenus=[
    go.layout.Updatemenu(
        buttons=buttons,
        direction="down",
        xanchor='left',
        yanchor='bottom'
    )],
    margin=go.layout.Margin(
        l=50,
        r=10,
        pad=4
    )
)
fig.show()
"""
