import plotly.graph_objects as go
import pandas as pd
import os

class plotLevel:
    @staticmethod
    def prepareDataframeDict(maindf):
        '''
        :param maindf: Dataframe of metagenome
        :param mapdf: Mapping dataframe
        :param metadavar: Variable from mapping dataframe
        :param levelName: Annotation level name or list of levels
        :param num: number of samples
        :return: dictionary of melted dataframes with keys as level names and values as melted dataframes
        '''
        df_dict = {}
        for l in set(maindf['level']):
            temp_df = maindf.loc[maindf['level'] == l]
            temp_df = temp_df.groupby(['taxon','level','sample']).sum().reset_index()
            temp_df = temp_df.drop(['level','taxon_abund','taxon_rel_abun','level_rel_abund','level_abund'],axis=1)
            temp_df = temp_df.pivot(index='taxon',columns='sample',values='taxon_contribution')
            temp_df = temp_df.div(temp_df.sum())*100
            temp_df = temp_df.loc[temp_df.sum(axis=1).sort_values(ascending=False)[0:10].index]
            temp_df.loc['Others'] = list(100-temp_df.sum())
            temp_df['taxon'] = temp_df.index
            temp_df = temp_df.melt(id_vars='taxon',value_name='taxon_contribution')
            temp_df = temp_df.sort_values(by='taxon_contribution',ascending=False)
            df_dict[l] = temp_df
        return df_dict
    @staticmethod
    def addTrace(df_dict):
        '''
        :param df_dict: Dictionary of dataframes to be plotted
        :return: Plotly.graph_object with traces
        '''
        fig = go.Figure()
        for k in sorted(df_dict.keys()):
            v = df_dict[k]
            df_temp = df_dict[k]
            nameList = list(set(df_temp['taxon']))
            for i, j in enumerate(nameList):
                x = df_temp[df_temp['taxon'] == j]['sample'].to_list()
                y = df_temp[df_temp['taxon'] == j]['taxon_contribution'].to_list()
                name = nameList[i]
                fig.add_trace(go.Bar(x=x, y=y, name=name, visible=False))
        fig.update_layout(barmode='stack')
        return fig
    @staticmethod
    def prepareAnnotations(df_dict):
        '''
        :param df_dict:  Dictionary of dataframes to be plotted
        :return:  Dictionary of annotations
        '''
        annotations = dict()
        for k in sorted(df_dict.keys()):
            df_temp = df_dict[k]
            temp = []
            for i, j in enumerate(set(df_temp['taxon'])):
                temp.append(
                    dict(
                        x=df_temp[df_temp['taxon'] == j]['sample'].to_list(),
                        y=df_temp[df_temp['taxon'] == j]['taxon_contribution'].to_list(),
                        text="Taxon contribution"
                    )
                )
            annotations[k] = temp
        return annotations
    @staticmethod
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
                    {'title': 'Taxon contribution for %s' % j,
                   'annotations': annotations[j]}]
                )
            buttons.append(temp)
            start = end
        return buttons
    @staticmethod
    def addUpdateMenus(fig, annotations, buttons):
        fig.update_layout(
            updatemenus=[
                go.layout.Updatemenu(
                    buttons=buttons,
                    direction="down",
                    xanchor='right',
                    yanchor='bottom'
                )],
            margin=go.layout.Margin(
                l=50,
                r=10,
                pad=4
            )
        )
        return fig

def plotDescription(df_file):
    df = pd.read_csv(df_file, sep='\t', compression='gzip')  # read input dataframe
    df['level'] = df['level'].str.replace('C\s\d+\s','')
    df['level'] = df['level'].str.replace('\[.*\]','')
    df['level'] = df['level'].str.strip()
    df = df[~df['level'].str.contains('Others')]
    # prepare dataframe (data to plot)
    df_dict = plotLevel.prepareDataframeDict(df)
    df_dict = {k:df_dict[k] for k in sorted(df_dict)}
    fig = plotLevel.addTrace(df_dict)
    annotations = plotLevel.prepareAnnotations(df_dict)
    buttons = plotLevel.prepareButtons(annotations,list(df_dict.keys()))
    fig = plotLevel.addUpdateMenus(fig, annotations, buttons)
    return fig