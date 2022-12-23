import os
import subprocess
import argparse

"""
Created on August 2021
@email: liaoshangfeng@genomics.cn
"""

parser = argparse.ArgumentParser(
    prog='SCAT Trajectory', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('Trajectory', help='Perform trajectory analysis by monocle3')

parser.add_argument(
    '-r', '--Rscript',
    dest='Rscript',
    type=str,
    default="/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/R411/bin/Rscript",
    help='Rscript path (default: %(default)s)')

parser.add_argument(
    '-i', '--input',
    dest='input',
    type=str,
    required=True,
    default=None,
    help='flt_norm_cluster_SeuratObject RDS file')

parser.add_argument(
    '-o', '--out',
    dest='out',
    type=str,
    required=True,
    default=None,
    help='Directory to save file')

parser.add_argument(
    '--prefix',
    type=str,
    default='sample',
    help='Sample ID, will be used as output prefix and seuratObj ident')

parser.add_argument(
    '-g', '--gene',
    dest='gene',
    type=str,
    default=None,
    help='gene name to define root of trajectory (default: %(default)s)')

parser.add_argument(
    '--mDEG',
    action='store_true',
    default=False,
    help='Add "--mDEG" to allow differential expression analysis (default: %(default)s)')

parser.add_argument(
    '--mGCN',
    action='store_true',
    default=False,
    help='Add "--mGCN" to allow gene co-expression analysis (default: %(default)s)')

parser.add_argument(
    '--mTop_markers',
    action='store_true',
    default=False,
    help='Add "--mTop_markers" to allow top markers analysis (default: %(default)s)')

args = parser.parse_args()


def get_r_script():
    r_script = os.path.join(os.path.dirname(__file__), "R/trajectory_monocle3_v2.R")
    return r_script


def make_dir(path):
    if(os.path.exists(path) and os.path.isdir(path)):
        pass
    else:
        os.makedirs(path)


def main():
    make_dir(args.out)
    r_script = get_r_script()
    commands = [args.Rscript, r_script]
    rscript_args = [
        "-i", args.input,
        "-o", args.out, 
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
    main()
