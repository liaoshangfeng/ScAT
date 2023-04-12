import os
import subprocess
import argparse

"""
Created on August 2021
@email: liaoshangfeng@genomics.cn
"""

def get_r_script():
    r_script = os.path.join(os.path.dirname(__file__), "R/qc_seurat_v10.R")
    return r_script


def make_dir(path):
    if(os.path.exists(path) and os.path.isdir(path)):
        pass
    else:
        os.makedirs(path)
        

def QcControlParser(subparsers):
    workflow = subparsers.add_parser("QcControl", help = "Perform QC control by Seurat")
    group_input = workflow.add_argument_group("input arguments")
    group_input.add_argument('--Rscript',type=str,required=True,default=None,help='Rscript path (default: %(default)s)')
    group_input.add_argument('--input', dest='input', type=str, required=True, default=None, help='Input Stereo-seq *.gem; 10X single cell directory; raw_matrix.txt.gz; raw_seuratObj.rds filename')
    group_input.add_argument('--binsize', dest='binsize', type=int, default=1, help='Bin size to binning, your input should be in binSize1 if you set this (default %(default)s)')
    group_input.add_argument('--geneColumn', dest='geneColumn', type=int, default=2, help='Specify which column of genes.tsv or features.tsv to use for gene names (default %(default)s)')
    group_input.add_argument('--fileType',dest='fileType',type=str,default='gem', choices=('gem', 'matrix'), help='Set as matrix to support *.matrix.txt.gz, otherwise, no need to change (default %(default)s)')

    group_output = workflow.add_argument_group("output arguments")
    group_output.add_argument('--out',type=str,required=True,default=None,help='Directory to save file')
    group_output.add_argument('--prefix', dest='prefix', type=str, default='sample', help='sample ID, will be used as output prefix and seurat object ident (default %(default)s)')
    group_output.add_argument('--ptsize', dest='ptsize', type=float, default=1, help='Point size for spatial dim plot (default %(default)s)')
    
    group_qc = workflow.add_argument_group("quality control arguments")
    group_qc.add_argument('--mtPattern', dest='mtPattern', type=str, default='^MT-', help='Mitochondria gene pattern (default %(default)s)')
    group_qc.add_argument('--rbPattern', dest='rbPattern', type=str, default='^RP[SL]', help='Ribosome gene pattern (default %(default)s)')
    group_qc.add_argument('--hbPattern', dest='hbPattern', type=str, default='^HB[^(P)]', help='Hemoglobin gene pattern (default %(default)s)')
    group_qc.add_argument('--platPattern', dest='platPattern', type=str, default='PECAM1|PF4', help='platelet gene pattern (default %(default)s)')
    
    group_filter = workflow.add_argument_group("filter arguements")
    group_filter.add_argument('--mFilter', dest='mFilter', action='store_true',default=False, help='Add "--mTop_markers" to allow filteration (default: %(default)s)')
    group_filter.add_argument('--minCount', dest='minCount', type=int, default=0, help='Minimum UMI number per cell/BIN [default %(default)s]')
    group_filter.add_argument('--maxCount', dest='maxCount', type=int, default=1000000, help='Maximum UMI number per cell/BIN [default %(default)s]')
    group_filter.add_argument('--minFeature', dest='minFeature', type=int, default=0, help='Minimum feature number per cell/BIN [default %(default)s]')
    group_filter.add_argument('--maxFeature', dest='maxFeature', type=int, default=100000, help='Maximum feature number per cell/BIN [default %(default)s]')
    group_filter.add_argument('--minCell', dest='minCell', type=int, default=3, help='Minimum cell or BIN number express the gene (for genes filterationn) (default %(default)s)') 
    group_filter.add_argument('--mtRatio', dest='mtRatio', type=float, default=100, help='Maximum mitochondria gene UMI ratio per cell/BIN (default %(default)s)')
    group_filter.add_argument('--rbRatio', dest='rbRatio', type=float, default=100, help='Maximum ribosome gene UMI ratio per cell/BIN (default %(default)s)')
    group_filter.add_argument('--hbRatio', dest='hbRatio', type=float, default=100, help='Maximum hemoglobin gene UMI ratio per cell/BIN (default %(default)s)')
    group_filter.add_argument('--platRatio', dest='platRatio',type=float, default=100, help='Maximum platelet gene UMI ratio per cell/BIN (default %(default)s)')


def QcControl(args):
    make_dir(args.out)
    r_script = get_r_script()
    commands = [args.Rscript, r_script]
    rscript_args = [
        "--input", args.input,
        "--out", args.out,
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
    parser = argparse.ArgumentParser(prog='ScAT', description='ScAT - a python xxxxxxxxxxxxxxxxxx',
                                  epilog='Use %(prog)s {command} -h to get help on individual commands')

    subparsers = parser.add_subparsers(dest = "subcommand")
    QcControlParser(subparsers)
    args = parser.parse_args()
    # print(args)
    QcControl(args)
