import os
from  math import floor,ceil
import argparse


parser = argparse.ArgumentParser(
    prog='SCAT TFregulon', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('TFregulon', help='Perform single-cell regulatory network inference')
parser.add_argument('--GeneMatrix',
        dest='genematrix',
        type=str,
        required=False,
        help='Gene expression matrix of single cell')

parser.add_argument('--species',
        dest='species',
        type=str,
        required=False,
        help='human or mouse')

parser.add_argument('--annot',
        dest='annot',
        type=str,
        required=False,
        help='Cell type annotation file')

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



def main():
    

    if os.path.exists(args.outdir)==False:
        os.system('mkdir '+args.outdir)
    
    ##调用human 脚本
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
                  ' --TFDataset '+os.path.join(args.tfdataset,args.species) +' --cpu '+ \
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
    main()



