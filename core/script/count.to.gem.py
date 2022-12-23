import pandas as pd

inputfile='ST.ms.brain.tsv'

outfile='./ST.ms.gem.tsv'

def matrix_gem():
    data=pd.read_table('ST.ms.brain.tsv',dtype=str)
    df_list=list()
    for col in data.columns:
        df=pd.DataFrame(data[col])
        df=df.loc[df[col]!='0']
        df=pd.DataFrame({'geneID':list(df.index),'x':col.split('_')[0],'y':col.split('_')[1],'MIDCounts':list(df[col])})
        df_list.append(df)
    df_concat=pd.concat(df_list)
    df_concat.to_csv('./ST.ms.gem.tsv',index=None,sep='\t')
matrix_gem()
