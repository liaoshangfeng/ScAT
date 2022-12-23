import os
from  math import floor,ceil
import argparse
import operator as op
from cytoolz import compose
import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import anndata as ad
from adjustText import adjust_text
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')
from pyscenic.utils import load_motifs
from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss
from pyscenic.utils import modules_from_adjacencies
import loompy as lp
from pyscenic.cli.utils import load_signatures
from matplotlib.backends.backend_pdf import PdfPages


parser = argparse.ArgumentParser()
parser.add_argument('--GeneMatrix',
        dest='genematrix',
        type=str,
        required=False,
        help='Gene expression matrix of single cell')

parser.add_argument('--annot',
        dest='annot',
        type=str,
        required=False,
        help='Cell type annotation file')


parser.add_argument('--species',
        dest='species',
        type=str,
        required=False,
        help='human or mouse')


parser.add_argument('--TFDataset',
        dest='tfdataset',
        type=str,
        required=False,
        help='TF dataset')
parser.add_argument('--TFName',
       dest='tf',
       type=str,
       required=False,
       help='Provide a TF name')
parser.add_argument('--NetFile',
        dest='net',
        type=str,
        required=False,
        help='Co-expression network inferred by pyscenic named \'02.Network.Adj.tsv\' in the outdir')

parser.add_argument('--SCENICloom',
        dest='ScenicLoom',
        type=str,
        required=False,
        help='Loom file outputed by SCENIC named \'SCENIC.Regulon.AUC.loom \'')



parser.add_argument('--h5ad',
        dest='h5ad',
        type=str,
        required=False,
        help='h5ad file including AdataExpression and regulonAUC value info named \'Adata.Expression.Regulon.AUC.h5ad\'')


parser.add_argument('--TFreglist',
       dest='TFreglist',
       type=str,
       required=False,
       help='TF regulons list file;One column')


parser.add_argument('--outdir',
        dest='outdir',
        type=str,
        required=True,
        help='Output result directory')

parser.add_argument('--cpu',
        dest='cpu',
        type=str,
        required=False,
        help='Cpu numbers')
args = parser.parse_args()

'''
01 函数功能:整合细胞类型文件和表达矩阵文件,产生loom文件
'''
def Generate_Expression_loom():
    df_annot = pd.read_table(args.annot)
    df_annot.columns=['CellID','CellType']
    df_tpm = pd.read_table(args.genematrix, index_col=0)
    adata = sc.AnnData(X=df_tpm.T.sort_index())
    df_obs = df_annot[['CellID','CellType']].set_index('CellID').sort_index()
    adata.obs = df_obs
    adata.var_names_make_unique()
    ##运行过慢，识别是否是测试数据；如果是，过滤大部分数据，加快测试数据后续计算
    ##实际数据不此函数不进行过滤
    if args.annot.split('/')[-1]=='GSE115978_cell.annotations.tsv':
        sc.pp.filter_cells(adata, min_genes=3500)
        sc.pp.filter_genes(adata, min_cells=3200)
    
    else:
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_cells(adata, min_genes=3)
    adata.raw = adata 
    #sc.pp.log1p(adata)
   
    row_attrs = { 
    "Gene": np.array(adata.var.index) ,
     }    
    col_attrs = { 
    
    "CellID":  np.array(adata.obs.index) ,
    'CellType':np.array(adata.obs['CellType']),
    'n_genes': np.array(adata.obs['n_genes']) ,
    }

    lp.create( os.path.join(args.outdir,'01.Adata.Expression.loom') ,adata.X.transpose(), row_attrs, col_attrs ) 


'''
02 函数功能:根据上一个函数产生的loom文件及转录因子注释文件，推断基因共表达文件
'''
def Net_Infer():
   order='pyscenic grn --num_workers '+ args.cpu + ' --output '+ os.path.join(args.outdir,'02.Network.Adj.tsv') + ' --method grnboost2 '+ os.path.join(args.outdir,'01.Adata.Expression.loom') +' '+ os.path.join(args.tfdataset,'hs_hgnc_tfs.txt')
   
   print('调用pyscenic grn命令'+'\t'+order)
   os.system(order)

'''
03 函数功能:根据转录起始位点信息，优化调控单元
参数:1,adj网络文件;2,人类tss注释文件；3，motif注释文件;4,表达文件;5，output;6,CPU
'''

def Candidate_Regulon():
    order='pyscenic ctx '+os.path.join(args.outdir,'02.Network.Adj.tsv')+' '+os.path.join(args.tfdataset,'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather ')+ ' \
        --annotations_fname  '+ os.path.join(args.tfdataset,'motifs-v9-nr.hgnc-m0.001-o0.0.tbl') +' \
         --expression_mtx_fname '+os.path.join(args.outdir,'01.Adata.Expression.loom') +' \
           --mode "dask_multiprocessing" \
           --output '+os.path.join(args.outdir,'03.Regulon.Refine.csv')+'  \
           --num_workers '+args.cpu+'  \
            --mask_dropouts'
    
    print('调用pyscenic ctx命令'+'\t'+order)
    os.system(order)

'''
04 函数功能:计算每个调控单元在每个细胞中的AUC值
##参数:1,表达文件;2,regulon文件;3,out;4,CPU
'''
def AUC_Regulon_Activity():
    order='pyscenic aucell '+os.path.join(args.outdir,'01.Adata.Expression.loom') +' '+ os.path.join(args.outdir,'03.Regulon.Refine.csv') +' '+'  \
         --output '+os.path.join(args.outdir,'04.SCENIC.Regulon.AUC.loom')+' \
          --num_workers '+args.cpu
    print('调用pyscenic aucell命令'+'\t'+order)
    os.system(order)

'''
05 函数功能:合并表达文件、调控单元文件和每个调控单元的AUC值
'''
def Merge_Expression_Regulon_AUC():
    adata=sc.read(os.path.join(args.outdir,'01.Adata.Expression.loom'))
    lf = lp.connect( os.path.join(args.outdir,'04.SCENIC.Regulon.AUC.loom'), mode='r', validate=False )
    auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
    sig = load_signatures(os.path.join(args.outdir,'03.Regulon.Refine.csv'))
    adata = add_scenic_metadata(adata, auc_mtx, sig)
    adata.write_h5ad(os.path.join(args.outdir,'05.Adata.Expression.Regulon.AUC.h5ad'))
    rss = regulon_specificity_scores(auc_mtx, adata.obs['CellType'])
    rss.to_csv(os.path.join(args.outdir,'06.Regulon.Cluster.RSS.value.tsv'),sep='\t')
    Plot_RSS_Group(rss)


'''
06 函数功能:定义的画图函数，Plot_RSS_Group函数调用
'''

def Rss_plot(rss, cell_type, top_n=5, max_n=None, ax=None):
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(4, 4))
    if max_n is None:
        max_n = rss.shape[1]
    data = rss.T[cell_type].sort_values(ascending=False)[0:max_n]
    #data=data.round(3)
    #ax=sns.scatter(
    ax.plot(np.arange(len(data)), data, '.',markersize=15)
    ax.set_ylim([(floor(data.min() * 100) / 100) ,(ceil(data.max() * 100) / 100)])
    ax.set_ylabel('RSS')
    ax.set_xlabel('Regulon')
    ax.set_title(cell_type)
    ax.set_xticklabels([])
    font = {'color': '#BA55D3', #'weight': 'normal',#'size': 4,
    }
    for idx, (regulon_name, rss_val) in enumerate(zip(data[0:top_n].index, data[0:top_n].values)):
        plt.style.use('classic')
        ax.plot([idx, idx], [rss_val, rss_val], 'r.')
        ax.text(
            idx + (max_n / 25),
            rss_val,
            regulon_name,
            fontdict=font,
            horizontalalignment='left',
            verticalalignment='center',
        )

'''
07 函数功能:画出每个cluster中top 5特异的调控单元
'''

def Plot_RSS_Group(rss):
    import matplotlib.ticker as ticker
    sns.set(style='whitegrid', font_scale=0.8)
    plt.rc('font',family='DejaVu Sans',size=16)
    #rss=pd.read_table(os.path.join(args.outdir,'05.Regulon.Cluster.RSS.value.tsv'), index_col=[0])
    cats = sorted( list(rss.index) )
    fig = plt.figure(figsize=(18, 12))
    #rss=rss.round(2)
    print('调控单元在每个cluster中的RSS值')
    print(rss)
    for c,num in zip(cats, range(1,len(cats)+1)):
        plt.style.use('classic')
        x=rss.T[c]
        ax = fig.add_subplot(2,ceil(len(rss)/2),num)
        Rss_plot(rss, c, top_n=5, max_n=None, ax=ax)
        ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
        for t in ax.texts:
            t.set_fontsize(17)
        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set_title(str(c),fontsize=20)
        ax.tick_params(axis='both', which='major', labelsize=14)
        #ax.get_yaxis().set_major_formatter(ticker.FormatStrFormatter('%.3f'))
        adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )
    fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='25')
    fig.text(0,0.5,'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='25')
    plt.tight_layout()
    #plt.rcParams.update({'figure.autolayout': True,})
    plt.savefig(os.path.join(args.outdir,"07.Cluster.RSS.top5.pdf"), bbox_inches = "tight")

'''
08 函数功能:默认输出2倍cluster数目的调控单元热图(为了做图好看)
'''
def Regulon_heatmap():
    plt.rcParams['figure.figsize'] = (15, 10)    

    pdf = PdfPages(os.path.join(args.outdir,'08.Cluster.RSS.heatmap.pdf'))
    rss_cluster=pd.read_table(os.path.join(args.outdir,'06.Regulon.Cluster.RSS.value.tsv'),index_col=[0])
    df_max = pd.DataFrame(rss_cluster.max(axis=0))
    ##按RSS最大值排序,默认挑选cluster 数目两倍的TF
    try:
        filter_TF=df_max.sort_values([0], ascending=False).iloc[len(rss_cluster)*2][0]
    except:
        filter_TF=df_max.sort_values([0], ascending=False).iloc[int(len(rss_cluster.columns)/2)][0]
    new_col = df_max.loc[df_max[0] > filter_TF].index
    rss_cluster=rss_cluster[new_col]
    col_list=rss_cluster.columns
    rss_cluster['cluster']=rss_cluster.index
    df_elem_list=list()
    for col in col_list:
       df_elem=rss_cluster[[col,'cluster']]
       df_elem.columns=['RSS_value','cluster']
       df_elem['regulon']=[col]*len(rss_cluster)
       df_elem_list.append(df_elem)
    cluster_rss=pd.concat(df_elem_list)
    df_heatmap = pd.pivot_table(data=cluster_rss,index='cluster', columns='regulon', values='RSS_value')
    fig, ax1 = plt.subplots()
    color=sns.color_palette("coolwarm", 10)
    sns.heatmap(df_heatmap.T,ax=ax1,annot=False, linewidths=.7, cbar=True, cbar_kws={"shrink": .45,'label':'RSS'},square=True, linecolor='gray',
                cmap=color, annot_kws={"size": 2})
    ax1.set_ylabel('')
    ax1.set_xlabel('')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    pdf.savefig(transparent=True)
    pdf.close()


'''
09 函数功能:指定一组调控单元,画此组调控单元在所有cluster中的特异值热图
'''
def Specify_Regulon_heatmap():
    plt.rcParams['figure.figsize'] = (15, 10)
    pdf = PdfPages(os.path.join(args.outdir,'09.Specify_Regulon.Cluster.RSS.heatmap.pdf'))
    rss_cluster=pd.read_table(os.path.join(args.outdir,'06.Regulon.Cluster.RSS.value.tsv'),index_col=[0])
    TFreg=pd.read_table(args.TFreglist)
    TFreg.columns=['TFreg']
    rss_cluster=rss_cluster[TFreg['TFreg']]
    #print('指定的一组调控单元在每个cluster中的特异值'+'\n'+rss_cluster)
    col_list=rss_cluster.columns
    rss_cluster['cluster']=rss_cluster.index
    df_elem_list=list()
    for col in col_list:
       df_elem=rss_cluster[[col,'cluster']]
       df_elem.columns=['RSS_value','cluster']
       df_elem['regulon']=[col]*len(rss_cluster)
       df_elem_list.append(df_elem)
    cluster_rss=pd.concat(df_elem_list)
    df_heatmap = pd.pivot_table(data=cluster_rss,index='cluster', columns='regulon', values='RSS_value')
    fig, ax1 = plt.subplots()
    color=sns.color_palette("coolwarm", 10)
    sns.heatmap(df_heatmap.T,ax=ax1,annot=False, linewidths=.7, cbar=True, cbar_kws={"shrink": .45,'label':'RSS'},square=True, linecolor='gray',
                cmap=color, annot_kws={"size": 3})
    ax1.set_ylabel('')
    ax1.set_xlabel('')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    pdf.savefig(transparent=True)
    pdf.close()

'''
10 函数功能: 指定一个TF，输出此TF调控网络图的原始文件
'''
def Certain_TF_regulon(tfname):
    from pyscenic.utils import modules_from_adjacencies
    res_dir=os.path.join(args.outdir,tfname)
    if os.path.exists(res_dir)==False:
        try:
            os.system('mkdir '+res_dir)
        except:
            pass

    adjacencies = pd.read_csv(args.net, index_col=False, sep='\t')
    lf = lp.connect(args.ScenicLoom, mode='r', validate=False)
    exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
    modules =modules_from_adjacencies(adjacencies, exprMat)
    modules=list(modules)

    regulons = {}
    for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
        regulons[i] =  list(r[r==1].index.values)

    tf = tfname
    tf_mods = [ x for x in modules if x.transcription_factor==tf ]

    for i,mod in enumerate( tf_mods ):
        print( f'{tf} module {str(i)}: {len(mod.genes)} genes' )
        print( f'{tf} regulon: {len(regulons[tf+"(+)"])} genes' )

    for i,mod in enumerate( tf_mods ):
        with open( os.path.join(res_dir,tf+'_module_'+str(i)+'.txt'), 'w') as f:
            for item in mod.genes:
                f.write("%s\n" % item)

    with open( os.path.join(res_dir,args.tf+'_regulon.txt'), 'w') as f:
        for item in regulons[tf+'(+)']:
            f.write("%s\n" % item)


'''
11 函数功能: 根据指定的TF，画图此TF在所有cluster中的热图
'''
def Plot_Certain_Regulon_AUC_umap(tfname):
    regulon_name='Regulon('+tfname+'(+))'
   # pdf= PdfPages('D:/BGI/04.pipeline/test.data.2/tsne.pdf')
    adata=sc.read_h5ad(args.h5ad)
    plt.rcParams['figure.figsize'] = (15, 10)
    sc.set_figure_params(frameon=False, fontsize=18)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=['CellType','Regulon(IRF1(+))'],use_raw=False,cmap='coolwarm', legend_fontsize=12) 
    
    plt.savefig(os.path.join(args.outdir,'10.'+regulon_name+'.Cluster.AUC.pdf'))


def main():
    import time
    import datetime
    if  os.path.exists(args.outdir)==False:
        try:
            os.system('mkdir '+args.outdir)
        except:
            pass

    start = datetime.datetime.now()
    now=int(round(time.time()*1000))
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(now/1000)))

    if args.tf:
        print('\n'+'开始获取'+args.tf+'的调控单元'+'\n')
        Certain_TF_regulon(args.tf)
        Plot_Certain_Regulon_AUC_umap(args.tf)
        end = datetime.datetime.now()
        print('\n'+'耗时 '+'\n'+str((end-start).seconds/60)+'分钟')
    elif args.TFreglist:
        print('输出特定TF调控单元的热图'+'\n') 
        Specify_Regulon_heatmap()

    else:
        print('整合表达数据和注释数据为loom文件'+'\n')
        Generate_Expression_loom()

        print('根据表达数据,推断共表达模块'+'\n')
        print('以下为 调用函数pyscenic grn的运行信息'+'\n')

        Net_Infer()
        print('推断潜在调控单元(regulon)'+'\n')
        Candidate_Regulon()
        print('计算每个调控单元在单个细胞中的AUC值;regulon在单个细胞中的活跃度'+'\n')
        AUC_Regulon_Activity()
        print('整合Cluster,AUC,regulon,GeneExpression信息'+'\n')
        Merge_Expression_Regulon_AUC()
        print('计算调控单元(regulon)在每个cluster中的RSS(Regulon specificity scores)'+'\n'+
              '范围(0,1),值越大代表在此cluster中越特异'+'\n')
        
        Regulon_heatmap()
        end = datetime.datetime.now()
        print('\n'+'耗时 '+'\n'+str(round((end-start).seconds/60,2))+'分钟')
if __name__ == "__main__":
    main()
    







