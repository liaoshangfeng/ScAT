import os
import argparse

parser = argparse.ArgumentParser(
    prog='SCAT DeconvSpot', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('DeconvSpot', help='Perform deconvolution of cell types by SPOTlight')

parser.add_argument('--Scdata',
        dest='scdata',
        type=str,
        required=True,
        help='Input seurate rds single-cell data after cell cluster annotation ')

parser.add_argument('--Stdata',
        dest='stdata',
        type=str,
        required=True,
        help='Input seurate rds spatial transcriptome data')

parser.add_argument('--piescale',
        dest='piescale',
        type=str,
        required=False,
        default=None,
        help='pie spot size ')

parser.add_argument('--ClusterFile',
        dest='clusterFile',
        type=str,
        required=False,
        default=None,
        help='View a set of specific cluster locations')

parser.add_argument('--Rscript',
        dest='rscript',
        type=str,
        required=False,
        default='/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/R411/bin/Rscript',
        help='The path of Rscript')

parser.add_argument('--outdir',
        dest='outdir',
        type=str,
        required=True,
        help='Output file into this directory')
args = parser.parse_args()

def r_code_filepath():
    filepath=os.path.join(os.path.split(os.path.realpath(__file__))[0],'R/DeconvSpot.r')
    return filepath

def main():
    if os.path.exists(args.outdir)==False:
        os.system('mkdir '+args.outdir)
    r_script=r_code_filepath()
    if args.clusterFile==None:
        order= args.rscript +' '+r_script +' --Scdata '+args.scdata +' --Stdata '+args.stdata +' --piescale '+args.piescale +' --outdir '+args.outdir
        print(order)
    else:
        order= args.rscript +' '+r_script +' --Scdata '+args.scdata +' --Stdata '+args.stdata +' --ClusterFile '+ args.clusterFile +' --piescale '+args.piescale +' --outdir '+args.outdir
        print(order)
    os.system(order)
    
if __name__ == "__main__":
    main()


