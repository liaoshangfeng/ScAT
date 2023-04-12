import sys
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
    r_script = os.path.join(os.path.dirname(__file__), "R/scCATCH_v1.R")
    return r_script


def AnnotscCATCHParser(subparsers):
    workflow = subparsers.add_parser("AnnotscCATCH", help = "Perform cell annotation by scCATCH")
    group_input = workflow.add_argument_group("input arguments")
    group_input.add_argument('--Rscript',type=str,required=True,default=None,help='Rscript path (default: %(default)s)')
    group_input.add_argument('--input', dest='input', type=str, required=True, default=None, help='Rds file of gene expression')
    group_input.add_argument('--species', dest='species', type=str, required=True, default=None,  help='Select species, support Human or Mouse (default: %(default)s)')
    group_input.add_argument('--tissue', dest='tissue', type=str, required=True, default=None,  help='tiisue type: https://github.com/ZJUFanLab/scCATCH/wiki (default: %(default)s)')
    group_input.add_argument('--cancer', dest='cancer', type=str, default=None,  help='cancer type:https://github.com/ZJUFanLab/scCATCH/wiki (default: %(default)s)')

    group_output = workflow.add_argument_group("output arguments")
    group_output.add_argument('--out',type=str,required=True,default=None,help='Directory to save file')
    group_output.add_argument('--prefix', dest='prefix', type=str, default='sample', help='sample ID, will be used as output file prefix (default: %(default)s)')


def AnnotscCATCH(args):
    if os.path.exists(args.out)==False:
        os.system('mkdir '+args.out)
    r_script = get_r_script()
    commands = [args.Rscript, r_script]
    if args.cancer==None:
        rscript_args = [
            "--input", args.input,
            '--tissue', args.tissue,
            '--species', args.species,
            "--out", args.out] 
    else:
        rscript_args = [
            "--input", args.input,
            "--cancer", args.cancer,
            "--tissue", args.tissue,
            "--species", args.species,
            "--out", args.out]    
            
    commands.extend(rscript_args)   
    subprocess.run([str(i) for i in commands])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ScAT', description='ScAT - a python xxxxxxxxxxxxxxxxxx',
                                  epilog='Use %(prog)s {command} -h to get help on individual commands')

    subparsers = parser.add_subparsers(dest = "subcommand")
    AnnotscCATCHParser(subparsers)
    args = parser.parse_args()
    AnnotscCATCH(args)


