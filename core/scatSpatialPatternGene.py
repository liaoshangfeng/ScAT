import os
import argparse

parser = argparse.ArgumentParser(
    prog='SCAT SpatialPatternGene', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('SpatialPatternGene', help='Statistical analysis of spatial expression pattern by SPARK')

parser.add_argument('--input',
        dest='input',
        type=str,
        required=True,
        help='Seurat rds after clustering  ')

parser.add_argument('--genefile',
        dest='genefile',
        type=str,
        required=False,
        help='One column file including genes you want to show')

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
       filepath=os.path.join(os.path.split(os.path.realpath(__file__))[0],'R/SpatialPatternGene.r')
       return filepath
def main():
    if os.path.exists(args.outdir)==False:
        os.system('mkdir '+args.outdir)
    r_script=r_code_filepath()
    if args.genefile==None:
      order= args.rscript+' '+r_script +' -i '+args.input + ' -o '+args.outdir
    else:
      order= args.rscript+' '+r_script +' -i '+args.input + ' -g ' +args.genefile + ' -o '+args.outdir
    print(order)
    os.system(order)
    
if __name__ == "__main__":
    main()





