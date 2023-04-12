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
    r_script = os.path.join(os.path.dirname(__file__), "R/SpatialCluster_v1.r")
    return r_script


def SpatialClusterParser(subparsers):
    workflow = subparsers.add_parser("SpatialCluster", help = "Perform spatial cluster by BayeSpaces")
    group_input = workflow.add_argument_group("input arguments")
    group_input.add_argument('--Rscript',type=str,required=True,default=None,help='Rscript path (default: %(default)s)')
    group_input.add_argument('--Stdata', dest='Stdata', type=str, required=True, default=None, help='Seurat rds before clustering...')
    group_input.add_argument('--clusterNumber', dest='clusterNumber', type=str, required=False, help='Set cluster number of BayesSpace spatial cluster')
    group_input.add_argument('--iterations', dest='iterations', type=str, required=True, help='iterations number of MCMC')
    
    group_output = workflow.add_argument_group("output arguments")
    group_output.add_argument('--out',type=str,required=True,default=None,help='Directory to save file')


def SpatialCluster(args):
    make_dir(args.out) 
    r_script = get_r_script()
    commands = [args.Rscript, r_script]
    rscript_args = [
        "--Stdata", args.Stdata,
        "--clusterNumber", args.clusterNumber,
        "--iterations", args.iterations,
        "--out", args.out] 

    commands.extend(rscript_args)
    subprocess.run([str(i) for i in commands])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ScAT', description='ScAT - a python xxxxxxxxxxxxxxxxxx',
                                  epilog='Use %(prog)s {command} -h to get help on individual commands')

    subparsers = parser.add_subparsers(dest = "subcommand")
    SpatialClusterParser(subparsers)
    args = parser.parse_args()
    SpatialCluster(args)
