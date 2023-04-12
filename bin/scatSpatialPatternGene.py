import os
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
    r_script = os.path.join(os.path.dirname(__file__), "R/SPARK.R")
    return r_script

        
def SpatialPatternGeneParser(subparsers):
    workflow = subparsers.add_parser("SpatialPatternGene", help = "Statistical analysis of spatial expression pattern by SPARK")
    group_input = workflow.add_argument_group("input arguments")     
    group_input.add_argument('--Rscript', type=str, required=True, default=None,help='Rscript path.(default: %(default)s)')
    group_input.add_argument('--input', type=str, required=True, dest='input', default=None, help='Seurat rds.')
    # group_input.add_argument('--genefile', type=str, dest='genefile', default=None, help='One column file including genes you want to show.')
    
    group_output = workflow.add_argument_group("output arguments")
    group_output.add_argument('--out', dest='out', type=str,required=True,default=None,help='Directory to save file')
    group_output.add_argument('--prefix', dest='prefix', type=str, default='sample', help='sample ID, will be used as output file prefix (default: %(default)s)')


def SpatialPatternGene(args):
    make_dir(args.out)
    r_script = get_r_script()
    commands = [args.Rscript, r_script]

    # if args.genefile==None:
    #     rscript_args = [
    #         "--input", args.input,
    #         "--out", args.out,
    #         "--prefix", args.prefix]
    # else:
    #     rscript_args = [
    #         "--input", args.input,
    #         "--genefile", args.genefile,
    #         "--out", args.out]
    
    rscript_args = [
            "--input", args.input,
            "--out", args.out,
            "--prefix", args.prefix]
    
    commands.extend(rscript_args)
    subprocess.run([str(i) for i in commands])
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ScAT', description='ScAT - a python xxxxxxxxxxxxxxxxxx',
                                  epilog='Use %(prog)s {command} -h to get help on individual commands')
    
    subparsers = parser.add_subparsers(dest = "subcommand")
    SpatialPatternGeneParser(subparsers)
    args = parser.parse_args()
    SpatialPatternGene(args)