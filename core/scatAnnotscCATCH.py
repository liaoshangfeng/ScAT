import os
import argparse

parser = argparse.ArgumentParser(
    prog='SCAT AnnotscCATCH', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('AnnotscCATCH', help='Perform cell annotaion by scCATCH')

parser.add_argument('--input',
        dest='input',
        type=str,
        required=True,
        help='rds file of gene expression ')

parser.add_argument('--species',
        dest='species',
        type=str,
        required=True,
        help='Human or Mouse')

parser.add_argument('--tissue',
        dest='tissue',
        type=str,
        required=True,
        help='tiisue type: https://github.com/ZJUFanLab/scCATCH/wiki')

parser.add_argument('--cancer',
        dest='cancer',
        type=str,
        required=False,
        help='cancer type:https://github.com/ZJUFanLab/scCATCH/wiki')

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
       filepath=os.path.join(os.path.split(os.path.realpath(__file__))[0],'R/scCATCH_v1.R')
       return filepath
def main():
    if os.path.exists(args.outdir)==False:
        os.system('mkdir '+args.outdir)
    r_script=r_code_filepath()
    if args.cancer==None:
      order= args.rscript+' '+r_script +' -i '+args.input + ' -t ' +args.tissue +' -s ' + args.species + ' -o '+args.outdir
    else:

      order= args.rscript+' '+r_script +' -i '+args.input + ' -t ' +args.tissue +' -s ' + args.species + ' -c ' +args.cancer + ' -o '+args.outdir
    print(order)
    os.system(order)
    
if __name__ == "__main__":
    main()





