import sys
import os
import subprocess
import argparse

"""
Created on August 2021
@email: liaoshangfeng@genomics.cn
"""

parser = argparse.ArgumentParser(
    prog='SCAT Annotation', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    'Annotation',
    help='Perform cell annnotation by SingleR')

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
    help='Input normalized data (*txt.gz)')

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
    '--ref_name',
    dest='ref_name',
    type=str,
    default='HumanPrimaryCellAtlasData',
    choices=[
        'BlueprintEncodeData',
        'DatabaseImmuneCellExpressionData',
        'HumanPrimaryCellAtlasData',
        'ImmGenData',
        'MonacoImmuneData',
        'MouseRNAseqData',
        'NovershternHematopoieticData'],
    help='Ref database name (default: %(default)s)')

parser.add_argument(
    '-c', '--cluster',
    dest='cluster',
    type=str,
    default=None,
    help='cluster file generated by Seurat workflow (default: %(default)s)')


args = parser.parse_args()


def get_r_script():
    r_script = os.path.join(os.path.dirname(__file__), "R/annotation_SingleR_v1.R")
    return r_script


def make_dir(path):
    if(os.path.exists(path) and os.path.isdir(path)):
        pass
    else:
        os.makedirs(path)


def get_ref_path(ref_name):
    database = [
        'BlueprintEncodeData',
        'DatabaseImmuneCellExpressionData',
        'HumanPrimaryCellAtlasData',
        'ImmGenData',
        'MonacoImmuneData',
        'MouseRNAseqData',
        'NovershternHematopoieticData']

    if ref_name in database:
        current_path = os.path.dirname(__file__)
        db_path = current_path.strip("core") + "docs/celldex/"
        ref = os.path.join(db_path, ref_name + '.rds')
        ref_col = os.path.join(db_path, ref_name + '_coldata.rds')
        return(ref, ref_col)
    else:
        parser.print_help(sys.stderr)
        print("ref_name should be selected from following: {database}\n".format(database="\t".join(database)))
        sys.exit(1)


def main():
    ref, ref_col = get_ref_path(args.ref_name)
    make_dir(args.out)
    r_script = get_r_script()
    commands = [args.Rscript, r_script]
    rscript_args = [
        "-i", args.input,
        "-o", args.out,
        "--prefix", args.prefix,
        "--ref", ref,
        "--ref_col", ref_col]

    if args.cluster:
        rscript_args.extend(["--cluster", args.cluster])

    commands.extend(rscript_args)
    subprocess.run([str(i) for i in commands])


if __name__ == "__main__":
    main()
