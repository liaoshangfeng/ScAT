import os
import subprocess
import argparse

"""
Created on August 2021
@email: liaoshangfeng@genomics.cn
"""

def get_r_script():
    r_script = os.path.join(os.path.dirname(__file__), "R/DeconvSpot.R")
    return r_script


def make_dir(path):
    if(os.path.exists(path) and os.path.isdir(path)):
        pass
    else:
        os.makedirs(path)


def DeconvSpotParser(subparsers):
    workflow = subparsers.add_parser("DeconvSpot", help = "Perform deconvolution of cell types by SPOTlight")
    group_input = workflow.add_argument_group("input arguments")
    group_input.add_argument('--Rscript', type=str, required=True,default=None,help='Rscript path (default: %(default)s)')
    group_input.add_argument('--sc_rds', dest='sc_rds', type=str, required=True, default=None, help='Input seurate rds single-cell data after cell cluster annotation.')
    group_input.add_argument('--st_rds', dest='st_rds', type=str, required=True, default=None, help='Input seurate rds spatial transcriptome data.')
    group_input.add_argument('--piescale', dest='piescale', type=str, default=None, help='Pie spot size.')
    group_input.add_argument('--ClusterFile', dest='clusterFile', type=str, default=None, help='View a set of specific cluster locations.')

    group_output = workflow.add_argument_group("output arguments")
    group_output.add_argument('--out', dest='out', type=str,required=True,default=None,help='Directory to save file')


def DeconvSpot(args):
    make_dir(args.out)
    r_script = get_r_script()
    commands = [args.Rscript, r_script]
    if args.clusterFile==None:
        rscript_args = [
            "--sc_rds", args.sc_rds,
            "--st_rds", args.st_rds,
            "--piescale", args.piescale,
            "--out", args.out]
    else:
        rscript_args = [
            "--sc_rds", args.sc_rds,
            "--st_rds", args.st_rds,
            "--ClusterFile", args.ClusterFile,
            "--piescale", args.piescale,
            "--out", args.out]
        
    commands.extend(rscript_args)
    subprocess.run([str(i) for i in commands])

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ScAT', description='ScAT - a python xxxxxxxxxxxxxxxxxx',
                                  epilog='Use %(prog)s {command} -h to get help on individual commands')
    
    subparsers = parser.add_subparsers(dest = "subcommand")
    DeconvSpotParser(subparsers)
    args = parser.parse_args()
    DeconvSpot(args)