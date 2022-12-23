import os
import argparse

parser = argparse.ArgumentParser(
    prog='SCAT SpatialCluster', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('SpatialCluster', help='Perform spatial cluster by BayeSpaces')

parser.add_argument('--Stdata',
        dest='Stdata',
        type=str,
        required=True,
        help='Seurat rds after clustering  ')

parser.add_argument('--clusterNumber',
        dest='clusterNumber',
        type=str,
        required=False,
        help='Set cluster number of BayesSpace spatial cluster')

parser.add_argument('--iterations',
        dest='iterations',
        type=str,
        required=True,
        help='iterations number of MCMC')

parser.add_argument('--outdir',
        dest='outdir',
        type=str,
        required=True,
        help='directory to save result file')

parser.add_argument('--Rscript',
        dest='rscript',
        type=str,
        required=False,
        default='/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/R411/bin/Rscript',
        help='The path of Rscript')


args = parser.parse_args()
def r_code_filepath():
       filepath=os.path.join(os.path.split(os.path.realpath(__file__))[0],'R/SpatialCluster_v1.r')
       return filepath
def main():
    if os.path.exists(args.outdir)==False:
        os.system('mkdir '+args.outdir)
    r_script=r_code_filepath()
    order= args.rscript+' '+r_script +' --Stdata '+args.Stdata + ' --clusterNumber ' +args.clusterNumber + ' --iterations '+ args.iterations +' --outdir '+ args.outdir

    print(order)
    os.system(order)
    
if __name__ == "__main__":
    main()





