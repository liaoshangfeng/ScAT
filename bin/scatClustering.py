import os
import subprocess
import argparse

"""
Created on August 2021
@email: liaoshangfeng@genomics.cn
"""


def get_r_script():
    r_script = os.path.join(os.path.dirname(__file__), "R/cluster_seurat.R")
    return r_script


def make_dir(path):
    if(os.path.exists(path) and os.path.isdir(path)):
        pass
    else:
        os.makedirs(path)


def ClusteringParser(subparsers):
    workflow = subparsers.add_parser("Clustering", help = "Perform cell cluster by Seurat")

    group_input = workflow.add_argument_group("input arguments")
    group_input.add_argument('--Rscript',type=str,required=True,default=None,help='Rscript path (default: %(default)s)')
    group_input.add_argument('--input', dest='input', type=str, required=True, default=None, help='Input *_seuratObj.rds with SCT assay')
    group_input.add_argument('--pc', type=int, default=30, help='PCA dimensionality for clustering. (default: %(default)s)')
    group_input.add_argument('--resolution', type=float, default=0.8, help='Resolution for clustering, use a value above (below) 1.0 if you want obtain a larger (smaller) number of communities. (default: %(default)s)')
    group_input.add_argument('--ptsize', type=float, default=1, help='Point size for spatil dim plot')
    group_input.add_argument('--topN', type=int, default=3, choices=range(1, 10), help='Select top N makers in each cluster for further visualization (default: %(default)s)')

    group_output = workflow.add_argument_group("output arguments")
    group_output.add_argument('--out',type=str,required=True,default=None,help='Directory to save file')
    group_output.add_argument('--prefix', dest='prefix', type=str, default='sample', help='sample ID, will be used as output file prefix (default: %(default)s)')


    group_boolean = workflow.add_argument_group("boolean arguments")
    group_boolean.add_argument('--findSpatialVar', dest='findSpatialVar', action='store_true', default=False, help='Add "--findSpatialVar" to find spatial variable features by Seurat (default: %(default)s)')

def Clustering(args):
    make_dir(args.out) 
    r_script = get_r_script()
    commands = [args.Rscript, r_script]
    rscript_args = [
        "-i", args.input,
        "-o", args.out,
        "--prefix", args.prefix,
        "--ptsize", args.ptsize,
        "--pc", args.pc,
        "--resolution", args.resolution,
        "--topN", args.topN] 

    if args.findSpatialVar:
        rscript_args.extend(["--findSpatialVar"])

    commands.extend(rscript_args)
    subprocess.run([str(i) for i in commands])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ScAT', description='ScAT - a python xxxxxxxxxxxxxxxxxx',
                                  epilog='Use %(prog)s {command} -h to get help on individual commands')

    subparsers = parser.add_subparsers(dest = "subcommand")
    ClusteringParser(subparsers)
    args = parser.parse_args()
    Clustering(args)
