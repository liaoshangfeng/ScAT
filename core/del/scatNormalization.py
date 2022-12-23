import os
import subprocess
import argparse

"""
Created on August 2021
@email: liaoshangfeng@genomics.cn
"""

parser = argparse.ArgumentParser(
    prog='SCAT Normalization',
    description="Description: \n \n " + "Input file: \
            1. marker gene file\
            2. pehnotype file\n" + "The function of this script is to run a binomial Generalized Linear Mixed Model by eagle\n",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('Normalization', help='Perform normalization process')

parser.add_argument(
    '-r', '--Rscript',
    dest='Rscript',
    type=str,
    default="/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/R411/bin/Rscript",
    help='Rscript path')

parser.add_argument(
    '-i', '--input',
    dest='input',
    type=str,
    required=True,
    default=None,
    help='Input raw data file (*.txt.gz)')

parser.add_argument(
    '-m', '--norm_method',
    dest='norm_method',
    type=str,
    default='Scran',
    choices=['Scran', 'Linnorm'],
    help='Select normalization method (default: %(default)s)')

parser.add_argument(
    '-o', '--out',
    dest='out',
    type=str,
    required=True,
    default=None,
    help='Directory to save file')


args = parser.parse_args()


def get_r_script():
    r_script = os.path.join(os.path.dirname(__file__), "R/normalization.R")
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
        "-m", args.norm_method,
        "-o", args.out]

    commands.extend(rscript_args)
    print(args)
    subprocess.run([str(i) for i in commands])


if __name__ == "__main__":
    main()
