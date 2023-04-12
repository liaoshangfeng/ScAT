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
    r_script = os.path.join(os.path.dirname(__file__), "R/trajectory_monocle3_v2.R")
    return r_script


def TrajectoryParser(subparsers):
    workflow = subparsers.add_parser("Trajectory", help = "Perform trajectory analysis by monocle3")
    group_input = workflow.add_argument_group("input arguments")
    group_input.add_argument('--Rscript', type=str, required=True, default=None, help='Rscript path (default: %(default)s)')
    group_input.add_argument('--input', dest='input', type=str, required=True, default=None, help='Input RDS file with seurat object, e.g. *flt_norm_cluster_SeuratObject.')
    group_input.add_argument('--gene', dest='gene', type=str,  default=None, help='Gene name to define root of trajectory (default: %(default)s)')           
    group_input.add_argument('--mDEG', dest='mDEG', action='store_true', default=False, help='Add "--mTop_markers" to allow filteration (default: %(default)s)')
    group_input.add_argument('--mGCN', dest='mGCN', action='store_true', default=False, help='Add "--mGCN" to allow gene co-expression analysis (default: %(default)s)')         
    group_input.add_argument('--mTop_markers', dest='mTop_markers', action='store_true', default=False, help='Add "--mTop_markers" to allow top markers analysis (default: %(default)s)')                                                                               
                             
    group_output = workflow.add_argument_group("output arguments")
    group_output.add_argument('--out',type=str,required=True,default=None,help='Directory to save file')
    group_output.add_argument('--prefix', type=str, default='sample', help='Sample ID, will be used as output prefix and seuratObj ident')
                        

def Trajectory(args):
    make_dir(args.out)
    r_script = get_r_script()
    commands = [args.Rscript, r_script]
    rscript_args = [
        "--input", args.input,
        "--out", args.out, 
        "--prefix", args.prefix]
                             
    if args.gene:
        rscript_args.extend(["--gene", args.gene])

    if args.mDEG:
        rscript_args.extend(["--mDEG"])

    if args.mGCN:
        rscript_args.extend(["--mGCN"])

    if args.mTop_markers:
        rscript_args.extend(["--mTop_markers"])

    commands.extend(rscript_args)
    subprocess.run([str(i) for i in commands])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ScAT', description='ScAT - a python xxxxxxxxxxxxxxxxxx',
                                  epilog='Use %(prog)s {command} -h to get help on individual commands')

    subparsers = parser.add_subparsers(dest = "subcommand")
    TrajectoryParser(subparsers)
    args = parser.parse_args()
    Trajectory(args)
