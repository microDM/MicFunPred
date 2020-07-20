import pandas as pd
import math 

df = pd.read_csv('/BioVolume/1_PhDWork/micfunPred_development/df_all_ko.tsv.gz',sep='\t',index_col=0)

genus = df.index.to_series().str.split('_',expand=True)[0].to_list()

final_df = pd.DataFrame(index=set(genus),columns=list(range(1,101,1)))

for i,tax in enumerate(set(genus)): 
    print(i,tax)
    # subset dataframe and sort by increasing genome size
    df_temp = df.loc[df.index.str.contains(tax)] 
    # calculate number of core genes at each percent
    for i in range(1,101,1): 
        l = df_temp.shape[0] 
        n = (i/100)*l 
        n = math.ceil(n)
        df_temp2 = df_temp[df_temp.columns[df_temp.astype(bool).sum()>=n]]
        final_df[i].loc[tax] = df_temp2.shape[1]
         
final_df.to_csv('geneCoverage.txt',sep='\t')
