#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Program:
#       This is the one click pipline for ScAT (scRNA-seq pipline)
# History:
# 2021/04/26   Junpu Mei(Author)  Firstrelease
# update:
#       rewrited the Main function of ScAT in ScAT file
# Author: Junpu Mei, meijunpu@genomics.cn; Yifei Sheng, shengyifei@bgi.cn
# Maintained by Shangfeng Liao, liaoshangfeng@genomics.cn


import os
import sys
import re
import time
import subprocess
import argparse

"""
Created on August 2021
@email: liaoshangfeng@genomics.cn
"""

parser = argparse.ArgumentParser(
    prog='SCAT Oneclick',
    description="Description: \n \n " + "Input file: \
            1. marker gene file\
            2. altCount file\
            3. pehnotype file\n" + "The function of this script is to run a binomial Generalized Linear Mixed Model by eagle\n",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('Oneclick', help='Perform single cell analysis by one click')

parser.add_argument(
    '-i', '--input',
    dest='input',
    type=str,
    required=True,
    default=None,
    help='input raw data (*.txt.gz)')

parser.add_argument(
    '-o', '--out',
    dest='out',
    type=str,
    required=True,
    default=None,
    help='Directory to save file')

parser.add_argument(
    '-c', '--config',
    dest='config',
    type=str,
    required=True,
    default=None,
    help='config.txt file')

args = parser.parse_args()


def read_config(config_file):
    parameter_config = {}
    if(os.path.isfile(config_file)):
        print(
            "\n[%s] processing configuration file: %s\n" %
            (time.asctime, config_file))
    else:
        sys.stderr.write("Error: %s does not exists!\n" % config_file)
        sys.exit(-1)

    with open(config_file, 'r') as f_in:
        for configinfo in f_in.readlines():
            if(re.search('#|=', configinfo)):
                # print("\t", configinfo.strip())
                if(re.search('\\|$', configinfo)):
                    print()
                if(re.search('=', configinfo)):
                    infolist = re.split("=", configinfo)
                    parameter_config[infolist[0].strip()] = infolist[1].strip()

    return parameter_config


def make_dir(path):
    if(os.path.exists(path) and os.path.isdir(path)):
        pass
    else:
        os.makedirs(path)


def get_file(path, suffix):
    filelist = []
    for parent, dirnames, filenames in os.walk(path):
        for filename in filenames:
            if filename.lower().endswith(suffix):
                filelist.append(os.path.join(parent, filename))
        return filelist


def initial_module(output_dir, module_name):
    all_path = [
        output_dir +
        "/process/" +
        module_name,
        output_dir +
        "/shell/" +
        module_name,
        output_dir +
        "/list/" +
        module_name,
        output_dir +
        "/upload/" +
        module_name]
    for path in all_path:
        make_dir(path)


def save_file(output_dir, module_name):
    wkd = output_dir + "/process/" + module_name

    images = get_file(wkd, ('.tiff', '.tif', '.pdf'))
    for image in images:
        subprocess.run(
            ['cp', image, output_dir + "/upload/" + module_name])

    tables = get_file(wkd, ('.txt', '.csv', '.txt.gz'))
    for table in tables:
        subprocess.run(
            ['cp', table, output_dir + "/list/" + module_name])

    scripts = get_file(wkd, ('.sh', '.r', '.R'))
    for script in scripts:
        subprocess.run(
            ['cp', script, output_dir + "/shell/" + module_name])


def normalize_matrix(matrix_filtered, output_dir, parameter_config):
    normalize_parameter = parameter_config['Normalize']
    rscript_path = parameter_config['Rscript']
    criteria = re.split('\\s', normalize_parameter)

    if(normalize_parameter == 'None'):
        matrix_normalized = matrix_filtered
        print(
            "[ %s ] Normalization module does not execute!\n\n" %
            time.asctime())
    else:
        initial_module(output_dir, "normalization")
        # configure parameter for normalization module
        wkd = output_dir + "/process/normalization"
        matrix_norm_name = os.path.basename(
            matrix_filtered).strip(".txt.gz") + "_norm.txt.gz"
        matrix_normalized = os.path.join(wkd, matrix_norm_name)
        norm_method = criteria[0]
        current_path = os.path.dirname(__file__)
        r_script = os.path.join(current_path, "normalization.R")

        # excute the r_script in directory /process/normalization
        subprocess.run(['cp', r_script, wkd])
        r_script_new = os.path.join(wkd, "normalization.R")
        norm_sh = os.path.join(wkd, "norm.sh")
        with open(norm_sh, 'w') as f_out:
            f_out.write("#!/usr/bin/bash\n")
            f_out.write(
                "{Rscript} {r_script} {wkd} {m1} {m2} {method} \n".format(
                    Rscript=rscript_path,
                    r_script=r_script_new,
                    wkd=wkd,
                    m1=matrix_filtered,
                    m2=matrix_norm_name,
                    method=norm_method))
        subprocess.run(['bash', norm_sh])

        save_file(output_dir, "normalization")
        print("[ %s ] Normalization module done!\n\n" % time.asctime())

    return matrix_normalized


def recognize_cell_type(matrix_in, output_dir, parameter_config, cluster_file):
    rscript_path = parameter_config['Rscript']
    recognize_parameter = parameter_config['CellType']
    criteria = re.split('\\s', recognize_parameter)

    initial_module(output_dir, "cell_annotation")
    # configure parameter for normalization module
    wkd = output_dir + "/process/cell_annotation"
    project_name = os.path.basename(
        matrix_in).strip(".txt.gz")
    cell_type = project_name + "_CellType.txt"
    celltype_file = os.path.join(wkd, cell_type)
    current_path = os.path.dirname(__file__)

    if criteria[0] == "SingleR":
        print(
            "[ %s ] Recognizing cell type, the tool and the parameter are as bellow:\n" %
            time.asctime())
        print("\tTools: %s" % criteria[0])
        print("\tParametaers: %s\n" % " ".join(criteria[1:]))
        r_script = current_path + '/SingleR.R'
        db_path = current_path.strip("core") + "docs/celldex/"
        ref_db = db_path + criteria[1] + ".rds"
        ref_db_col = db_path + criteria[1] + "_coldata.rds"

        # excute the r_script in directory /process/cell_annotation
        subprocess.run(['cp', r_script, wkd])
        r_script_new = os.path.join(wkd, "SingleR.R")

        anno_sh = os.path.join(wkd, "anno.sh")
        with open(anno_sh, 'w') as f_out:
            f_out.write("#!/usr/bin/bash\n")
            f_out.write(
                "{Rscript} {r_script} {wkd} {m} {ref} {ref_col} {cell_type} {cluster_file} \n".format(
                    Rscript=rscript_path,
                    r_script=r_script_new,
                    wkd=wkd,
                    m=matrix_in,
                    ref=ref_db,
                    ref_col=ref_db_col,
                    cell_type=celltype_file,
                    cluster_file=cluster_file))
        subprocess.run(['bash', anno_sh])
        save_file(output_dir, "cell_annotation")
    print("[ %s ] Cell annotation module done!\n\n" % time.asctime())

    return celltype_file


def cell_cluster(matrix_imputed, output_dir, parameter_config):
    cluster_parameter = parameter_config['Cluster']
    rscript_path = parameter_config['Rscript']
    criteria = re.split('\\s', cluster_parameter)
    cluster_file = None
    markers_file = None

    if(cluster_parameter == 'None'):
        print(
            "[ %s ] Clustering module does not execute!\n\n" %
            time.asctime())
    else:
        print(
            "[ %s ] processing clustering module, the tool and the parameter are as bellow:\n" %
            time.asctime())
        print("\tTools: %s" % criteria[0])

        initial_module(output_dir, "cell_cluster")
        # configure parameter for cell cluster module
        wkd = output_dir + "/process/cell_cluster"
        cluster_resolution = criteria[1]
        pca_dim = criteria[2]
        project_name = os.path.basename(
            matrix_imputed).strip(".txt.gz")
        current_path = os.path.dirname(__file__)
        r_script = os.path.join(current_path, "cluster_by_seurat.R")

        # excute the r_script in directory /process/normalization
        subprocess.run(['cp', r_script, wkd])
        r_script_new = os.path.join(wkd, "cluster_by_seurat.R")
        cluster_sh = os.path.join(wkd, "cluster.sh")
        with open(cluster_sh, 'w') as f_out:
            f_out.write("#!/usr/bin/bash\n")
            f_out.write(
                "{Rscript} {r_script} {wkd} {data} {prefix} {res} {dim}\n".format(
                    Rscript=rscript_path,
                    r_script=r_script_new,
                    wkd=wkd,
                    data=matrix_imputed,
                    prefix=project_name,
                    res=cluster_resolution,
                    dim=pca_dim))

        subprocess.run(['bash', cluster_sh])

        cluster_file = os.path.join(wkd, project_name + "_clusterInfo.txt")
        print(cluster_file)
        markers_file = os.path.join(wkd, project_name + "_markerGene.txt")

        celltype_file = recognize_cell_type(
            matrix_imputed, output_dir, parameter_config, cluster_file)

        save_file(output_dir, "cell_cluster")
        print("[ %s ] Cell cluster module done!\n\n" % time.asctime())

    return (cluster_file, markers_file, celltype_file)


def cell_cell_interaction(
        matrix_filtered,
        celltype_file,
        output_dir,
        parameter_config):
    cci_parameter = parameter_config['CCI']
    rscript_path = parameter_config['Rscript']

    if(cci_parameter == 'None'):
        print(
            "[ %s ] CCIAnalysis module does not execute!\n\n" %
            time.asctime())
    else:
        criteria = re.split('\\s', cci_parameter)
        print(
            "[ %s ] processing CCIAnalysis module, the tool and the parameter are as bellow:\n" %
            time.asctime())
        print("\tTools: %s" % criteria[0])
        print("\tParametaers: %s\n" % " ".join(criteria[1:]))

        initial_module(output_dir, "CCI")
        # configure parameter for CCI module
        wkd = output_dir + "/process/CCI"
        project_name = os.path.basename(
            matrix_filtered).strip(".txt.gz")
        current_path = os.path.dirname(__file__)
        if criteria[0] == "CSOmap":
            r_script = current_path + '/CSOmap_arg.R'
            lr_pairs = current_path.strip(
                "core") + "docs/CSOmap.R/LR_pairs.txt"

            print(
                "[ %s ] Recognizing cell type, the tool and the parameter are as bellow:\n" %
                time.asctime())

            # excute the r_script in directory /process/normalization
            subprocess.run(['cp', r_script, wkd])
            r_script_new = os.path.join(wkd, "CSOmap_arg.R")
            cci_sh = os.path.join(wkd, "cci.sh")
            with open(cci_sh, 'w') as f_out:
                f_out.write("#!/usr/bin/bash\n")
                f_out.write(
                    "{Rscript} {r_script} {wkd} {m} {cell_type} {lr_pairs}\n".format(
                        Rscript=rscript_path,
                        r_script=r_script_new,
                        wkd=wkd,
                        m=matrix_filtered,
                        cell_type=celltype_file,
                        lr_pairs=lr_pairs))

            subprocess.run(['bash', cci_sh])

            count_file = os.path.join(wkd, "counts.txt")
            p_value_file = os.path.join(wkd, "pvalue.txt")

            r_script = current_path + '/CCI_Visualization.R'
            # excute the r_script in directory /process/normalization
            subprocess.run(['cp', r_script, wkd])
            r_script_new = os.path.join(wkd, "CCI_Visualization.R")
            plot_sh = os.path.join(wkd, "plot.sh")
            with open(plot_sh, 'w') as f_out:
                f_out.write("#!/usr/bin/bash\n")
                f_out.write(
                    "{Rscript} {r_script} {wkd} {count_file} {p_value_file} {project_name}\n".format(
                        Rscript=rscript_path,
                        r_script=r_script_new,
                        wkd=wkd,
                        count_file=count_file,
                        p_value_file=p_value_file,
                        project_name=project_name))
            subprocess.run(['bash', plot_sh])

            save_file(output_dir, "CCI")


def main():
    make_dir(args.out)
    initial_module(args.out, "")
    parameter_config = read_config(args.config)
    matrix_normalized = normalize_matrix(
        args.input, args.out, parameter_config)
    cluster_file, markers_file, celltype_file = cell_cluster(
        matrix_normalized, args.out, parameter_config)
    cell_cell_interaction(
        args.input,
        celltype_file,
        args.out,
        parameter_config)


if __name__ == "__main__":
    main()
