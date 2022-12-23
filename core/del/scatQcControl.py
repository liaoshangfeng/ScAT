import os
import subprocess
import argparse

"""
Created on August 2021
@email: liaoshangfeng@genomics.cn
"""

parser = argparse.ArgumentParser(
    prog='SCAT QcControl', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('QcControl', help='Perform QC control by Seurat')

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
    help='input Stereo-seq *.gem; 10X single cell directory; raw_matrix.txt.gz; raw_seuratObj.rds filename')

parser.add_argument(
    '-o', '--out',
    dest='out',
    type=str,
    required=True,
    default=None,
    help='Directory to save file')

parser.add_argument(
    '-b', '--binsize',
    dest='binsize',
    type=int,
    default=1,
    help='bin size to binning, your input should be in binSize1 if you set this (default %(default)s)')

parser.add_argument(
    '--geneColumn',
    dest='geneColumn',
    type=int,
    default=2,
    help='Specify which column of genes.tsv or features.tsv to use for gene names (default %(default)s)')

parser.add_argument(
    '--fileType',
    dest='fileType',
    type=str,
    default='gem',
    choices=('gem', 'matrix'),
    help='Set as matrix to support *.matrix.txt.gz, otherwise, no need to change (default %(default)s)')

parser.add_argument(
    '--ptsize',
    dest='ptsize',
    type=float,
    default=1,
    help='point size for spatial dim plot (default %(default)s)')

parser.add_argument(
    '--prefix',
    dest='prefix',
    type=str,
    default='sample',
    help='sample ID, will be used as output prefix and seurat object ident (default %(default)s)')

qc_metrics = parser.add_argument_group("QC metrics")
qc_metrics.add_argument(
    '--mtPattern',
    dest='mtPattern',
    type=str,
    default='^MT-',
    help='mitochondria gene pattern (default %(default)s)')

qc_metrics.add_argument(
    '--rbPattern',
    dest='rbPattern',
    type=str,
    default='^RP[SL]',
    help='ribosome gene pattern (default %(default)s)')

qc_metrics.add_argument(
    '--hbPattern',
    dest='hbPattern',
    type=str,
    default='^HB[^(P)]',
    help='hemoglobin gene pattern (default %(default)s)')

qc_metrics.add_argument(
    '--platPattern',
    dest='platPattern',
    type=str,
    default='PECAM1|PF4',
    help='platelet gene pattern (default %(default)s)')

filter_module = parser.add_argument_group("Filter module")
filter_module.add_argument(
    '--mFilter',
    dest='mFilter',
    action='store_true',
    default=False,
    help='Add "--mTop_markers" to allow filteration (default: %(default)s)')

filter_module.add_argument(
    '--minCount',
    dest='minCount',
    type=int,
    default=0,
    help='minimum UMI number per cell/BIN [default %(default)s]')

filter_module.add_argument(
    '--maxCount',
    dest='maxCount',
    type=int,
    default=1000000,
    help='maximum UMI number per cell/BIN [default %(default)s]')

filter_module.add_argument(
    '--minFeature',
    dest='minFeature',
    type=int,
    default=0,
    help='minimum feature number per cell/BIN [default %(default)s]')

filter_module.add_argument(
    '--maxFeature',
    dest='maxFeature',
    type=int,
    default=100000,
    help='maximum feature number per cell/BIN [default %(default)s]')

filter_module.add_argument(
    '--minCell',
    dest='minCell',
    type=int,
    default=3,
    help='minimum cell or BIN number express the gene (for genes filterationn) (default %(default)s)')

filter_module.add_argument(
    '--mtRatio',
    dest='mtRatio',
    type=float,
    default=100,
    help='maximum mitochondria gene UMI ratio per cell/BIN (default %(default)s)')

filter_module.add_argument(
    '--rbRatio',
    dest='rbRatio',
    type=float,
    default=100,
    help='maximum ribosome gene UMI ratio per cell/BIN (default %(default)s)')

filter_module.add_argument(
    '--hbRatio',
    dest='hbRatio',
    type=float,
    default=100,
    help='maximum hemoglobin gene UMI ratio per cell/BIN (default %(default)s)')

filter_module.add_argument(
    '--platRatio',
    dest='platRatio',
    type=float,
    default=100,
    help='maximum platelet gene UMI ratio per cell/BIN (default %(default)s)')

args = parser.parse_args()


def get_r_script():
    r_script = os.path.join(os.path.dirname(__file__), "R/qc_seurat_v7.R")
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
        "--binsize", args.binsize,
        "--geneColumn", args.geneColumn,
        "--fileType", args.fileType,
        "--mtPattern", args.mtPattern,
        "--rbPattern", args.rbPattern,
        "--hbPattern", args.hbPattern,
        "--platPattern", args.platPattern,
        "--minCount", args.minCount,
        "--maxCount", args.maxCount,
        "--minFeature", args.minFeature,
        "--maxFeature", args.maxFeature,
        "--minCell", args.minCell,
        "--mtRatio", args.mtRatio,
        "--rbRatio", args.rbRatio,
        "--hbRatio", args.hbRatio,
        "--platRatio", args.platRatio]

    if args.mFilter:
        rscript_args.extend(["--mFilter"])

    commands.extend(rscript_args)
    subprocess.run([str(i) for i in commands])


if __name__ == "__main__":
    main()
