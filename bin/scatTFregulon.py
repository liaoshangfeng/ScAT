import os
from  math import floor,ceil
import subprocess
import argparse


"""
Created on August 2021
@email: liaoshangfeng@genomics.cn
"""

def make_dir(path):
    if(os.path.exists(path) and os.path.isdir(path)):
        pass
    else:
        os.makedirs(path)


def get_r_script():
    r_script = os.path.join(os.path.dirname(__file__), "R/SpatialCluster_v1.r")
    return r_script


def TFregulonParser(subparsers):
    workflow = subparsers.add_parser("TFregulon", help = "Perform single-cell regulatory network inference.")
    group_input = workflow.add_argument_group("input arguments")
    group_input.add_argument('--GeneMatrix', dest='genematrix', type=str, required=True, default=None, help='Gene expression matrix of single cell.')
    group_input.add_argument('--species', dest='species', type=str, required=True, default=None, help='Human or mouse.')
    group_input.add_argument('--annot', dest='annot', type=str, required=False, default=None, help='Cell type annotation file.')
    group_input.add_argument('--TFDataset', dest='tfdataset', type=str, required=False, default=None, help='TF dataset.')
    group_input.add_argument('--TFName', dest='tf', type=str, required=False, default=None, help='Provide a TF name.')
    group_input.add_argument('--NetFile', dest='net', type=str, required=False, default=None, help='Co-expression network inferred by pyscenic named \'02.Network.Adj.tsv\' in the outdir.')
    group_input.add_argument('--SCENICloom', dest='ScenicLoom', type=str, required=False, default=None, help='Loom file outputed by SCENIC named \'SCENIC.Regulon.AUC.loom \'.')
    group_input.add_argument('--h5ad', dest='h5ad', type=str, required=False, default=None, help='h5ad file including AdataExpression and regulonAUC value info named \'Adata.Expression.Regulon.AUC.h5ad\'.')
    group_input.add_argument('--TFreglist', dest='TFreglist', type=str, required=False, default=None, help='TF regulons list file;One column.')
    group_input.add_argument('--cpu', dest='cpu', type=str, required=False, default=None, help='Cpu numbers.')
    
    group_output = workflow.add_argument_group("output arguments")
    group_output.add_argument('--outdir', type=str, required=True, default=None,help='Output result directory')
    

def TFregulon(args):
    make_dir(args.outdir) 
    r_script = get_r_script()
    commands = [args.Rscript, r_script]
    
    ## use human 
    if args.species=='human':
        filepath=os.path.join(os.path.split(os.path.realpath(__file__))[0],'pycode/TFregulon.hum.py')
        if args.tfdataset is not None: 
            order='python3 '+ filepath +' --species '+args.species +' --GeneMatrix '+ args.genematrix +' --annot '+ args.annot + \
            '       --TFDataset '+os.path.join(args.tfdataset,args.species) +' --cpu '+ \
                          args.cpu +' --outdir '+ args.outdir  
            print(order)
            os.system(order)
        if args.TFreglist is not None:
            order=' python3 '+ filepath +' --species '+ args.species +' --TFreglist '+args.TFreglist +' --outdir '+ args.outdir 
            print(order)
            os.system(order)

        if args.tf is not None:
            order=' python3 '+ filepath +' --species '+ args.species +' --TFName '+args.tf +' --SCENICloom '+ args.ScenicLoom  +' \
                    --h5ad  '+ args.h5ad +' --NetFile '+args.net + ' --outdir '+ args.outdir 
            print(order)
            os.system(order)

            
    ##调用mouse脚本
    elif args.species=='mouse':
        filepath=os.path.join(os.path.split(os.path.realpath(__file__))[0],'pycode/TFregulon.ms.py')
        if args.tfdataset is not None:
            order='python3 '+ filepath +' --GeneMatrix '+ args.genematrix +' --annot '+ args.annot + \
                  ' --TFDataset '+os.path.join(args.tfdataset,args.species) + ' --species '+args.species +' --cpu '+ \
                          args.cpu +' --outdir '+ args.outdir  
            print(order)
            os.system(order)
        if args.TFreglist is not None:
            order=' python3 '+ filepath +' --species '+ args.species +' --TFreglist '+args.TFreglist +' --outdir '+ args.outdir 
            print(order)
            os.system(order)

        if args.tf is not None:
            order=' python3 '+ filepath +' --species '+ args.species +' --TFName '+args.tf +' --SCENICloom '+ args.ScenicLoom  +' \
                    --h5ad  '+ args.h5ad +' --NetFile '+args.net + ' --outdir '+ args.outdir 
            print(order)
            os.system(order) 
    else:
         print('\n'+'Error: please choose species: human or mouse'+'\n') 
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ScAT', description='ScAT - a python xxxxxxxxxxxxxxxxxx',
                                  epilog='Use %(prog)s {command} -h to get help on individual commands')

    subparsers = parser.add_subparsers(dest = "subcommand")
    TFregulonParser(subparsers)
    args = parser.parse_args()
    TFregulon(args)