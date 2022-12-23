import os
import subprocess
import argparse

"""
Created on August 2021
@email: liaoshangfeng@genomics.cn
"""

parser = argparse.ArgumentParser(
    prog='SCAT Clustering', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('Clustering', help='Perform seurat workflow')

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
    help='Input 10x single cell data directory; raw_matrix.txt.gz; stereo *.gem; raw_seuratObj.rds')

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
    '--pc',
    type=int,
    default=30,
    help='PCA dimensionality for clustering. (default: %(default)s)')

parser.add_argument(
    '--resolution',
    type=float,
    default=0.8,
    help='Resolution for clustering, use a value above (below) 1.0 if you want obtain a larger (smaller) number of communities. (default: %(default)s)')

parser.add_argument(
    '--ptsize',
    type=float,
    default=1,
    help='Point size for spatil dim plot')

parser.add_argument(
    '--topN',
    type=int,
    default=3,
    choices=range(1, 10),
    help='Select top N makers in each cluster for further visualization (default: %(default)s)')

parser.add_argument(
    '--mSpatialVar',
    dest='mSpatialVar',
    action='store_true',
    default=False,
    help='Add "--mSpatialVar" to find spatial variable (default: %(default)s)')

args = parser.parse_args()


def get_r_script():
    r_script = os.path.join(os.path.dirname(__file__), "R/cluster_seurat.R")
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
        "--prefix", args.prefix,
        "--ptsize", args.ptsize,
        "--pc", args.pc,
        "--resolution", args.resolution,
        "--topN", args.topN] 

    if args.mSpatialVar:
        rscript_args.extend(["--mSpatialVar"])

    commands.extend(rscript_args)
    subprocess.run([str(i) for i in commands])


if __name__ == "__main__":
    main()
