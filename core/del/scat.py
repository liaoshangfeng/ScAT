#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Program:
#       This is the pipline for ScAT (scRNA-seq pipline)
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
import numpy
import subprocess
import gzip

# --  import customized modules  -- #

# from . import imputation as im
from . import trajectory as tr
from . import enrichment as er


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


def filter_matrix(matrix_raw, output_dir, parameter_config):
    matrix_filtered_name = os.path.basename(
        matrix_raw).strip(".txt.gz") + "_filtered.txt.gz"
    filter_parameter = parameter_config['Filter']

    if(filter_parameter == 'None'):
        matrix_filtered = matrix_raw
        print("[ %s ] Filter module does not execute!\n\n" % time.asctime())
    else:
        initial_module(output_dir, "filtration")
        # configure parameter for filtration module
        matrix_filtered = os.path.join(
            output_dir, "process", "filtration", matrix_filtered_name)
        criteria = re.split('\\s', filter_parameter)
        print(criteria)
        print(
            "[ %s ] processing filter module, and the filter criteria are as bellow:\n" %
            time.asctime())
        print("\tgene expressed: TPM > %s" % criteria[0])
        print("\teffective gene: expressed in at least %s cells" %
              (str(100 * float(criteria[1])) + '%'))
        print("\texpression variable coefficient: %s" % criteria[2])
        print("\tfolds of standard deviation: %s\n" % criteria[3])

        cellDict = dict()
        geneDict = dict()
        rownum = 0
        listnum = 0
        rare_gene_rm = 0
        stable_gene_rm = 0
        final_gene_num = 0

        # read the original matrix
        with gzip.open(matrix_raw, 'rt') as f_in:
            for original_matrix_info in f_in.readlines():
                if not original_matrix_info.strip() == '':
                    original_matrix_info_List = re.split(
                        r'(?:[,\t\s+])', original_matrix_info.strip())
                    tmp = original_matrix_info_List[:]
                    for i in tmp:
                        if i == '':
                            original_matrix_info_List.remove(i)
                    del tmp

                    rownum += 1

                    if(rownum == 1):
                        listnum = len(original_matrix_info_List) - 1
                        geneDict[rownum] = "\t".join(original_matrix_info_List)
                        # print (geneDict[rownum])
                    else:
                        isexp = 0
                        isstable = 1
                        # expcutoff = criteria[0] * listnum
                        expnum = 0
                        expsum = 0
                        for g in (range(1, listnum + 1)):
                            if(re.search(r'[^0-9.]', original_matrix_info_List[g])):
                                sys.stderr.write(
                                    "Error: all expression data must be digital!\n")
                                sys.exit(-1)
                            else:
                                if original_matrix_info_List[g] > criteria[0]:
                                    expnum += 1
                                expsum += float(original_matrix_info_List[g])
                        if(expnum > listnum * float(criteria[1])):
                            isexp = 1
                        else:
                            rare_gene_rm += 1

                        expaverage = expsum / listnum
                        expstd = numpy.std(
                            list(map(float, original_matrix_info_List[1:])))
                        cv = 10000
                        if expaverage != 0:
                            cv = expstd / expaverage
                        if(cv < float(criteria[2])):
                            isstable = 0
                        else:
                            stable_gene_rm += 1
                        isstable = 0
                        # exclude genes with low expression rate or with stable
                        # expression level in all cells.
                        if(isstable == 0 and isexp == 1):
                            geneDict[rownum] = "\t".join(
                                original_matrix_info_List)
                            final_gene_num += 1
                            for c in (range(1, listnum + 1)):
                                if(original_matrix_info_List[c] > criteria[0]):
                                    if(c in cellDict):
                                        cellDict[c] += 1
                                    else:
                                        cellDict[c] = 1

        # exclude cells with abnormal number of expressed gene
        cell_std = numpy.std(list(map(float, cellDict.values())))
        cell_ave = numpy.average(list(map(float, cellDict.values())))
        cell_min = cell_ave - float(criteria[3]) * cell_std
        cell_max = cell_ave + float(criteria[3]) * cell_std
        cell_rm = 0
        final_cell_num = 0
        for c in (range(1, listnum + 1)):
            if(c in cellDict and (cellDict[c] < cell_min or cellDict[c] > cell_max)):
                del cellDict[c]
                cell_rm += 1
            else:
                final_cell_num += 1

        final_gene_num = 0

        # write out the new matrix
        with gzip.open(matrix_filtered, 'wt') as f_out:
            for g in range(rownum):
                if g in geneDict:
                    final_gene_num += 1
                    exp_List = re.split(r'(?<!\t)\t|\t(?!\t)', geneDict[g])
                    new_info = exp_List[0]

                    # The number of elements in the header line and the
                    # expression line differs by 1
                    num = listnum - 1
                    if g != 1:
                        num = listnum

                    for c in range(1, num):
                        if(c in cellDict):
                            new_info = new_info + "\t" + exp_List[c]
                    new_info = new_info + "\n"
                    f_out.write(new_info)

        # save the norm.data to /list/filtration
        save_file(output_dir, "filtration")

        print(
            "\tThe original expression matrix have %s genes and %s cells." %
            (rownum - 1, listnum))
        print("\t%s genes were discarded for expressed (TPM >= %s) in less than %s cells." % (
            rare_gene_rm, criteria[0], (str(100 * float(criteria[1])) + '%')))
        # print("\t%s genes were discarded for having stable expression
        # level(variable coefficient < %s) in all cells." % (stable_gene_rm,
        # criteria[1]))
        print(
            "\t%s cells were discarded for having too much(> %s) or too little(< %s) expressed gene." %
            (cell_rm, int(cell_max), int(cell_min)))
        print(
            "\tThe new expression matrix have %s genes and %s cells, and is listed in %s.\n" %
            (final_gene_num - 1, final_cell_num - 1, matrix_filtered))
        print("[ %s ] Filter module done!\n\n" % time.asctime())

    return matrix_filtered


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


def impute_missing_value(matrix_normalized, output_dir, parameter_config):
    impute_parameter = parameter_config['Imputation']
    if(impute_parameter == 'None'):
        matrix_imputed = matrix_normalized
        print("[ %s ] Imputaion module does not execute!\n\n" % time.asctime())
    else:
        # if not os.path.exists(output_dir +
        # "/imputation"):os.makedirs(output_dir + "/imputation")
        print("Error, imputation module is developing. Please set as \"None\" ")
#         matrix_imputed_name = os.path.basename(matrix_normalized).strip(".txt") + "_imputed.txt"
#         matrix_imputed = os.path.join(output_dir, "imputation", matrix_imputed_name)
#         criteria = re.split('\\s', impute_parameter)
#         print(
#             "[ %s ] processing Imputation module, the tool and the parameter are as bellow:\n" %
#             time.asctime())
#         print("\tTools: %s" % criteria[0])
#         print("\tParametaers: %s\n" % " ".join(criteria[1:]))
#         script_path = os.path.join(output_dir, "shell")
#         if(criteria[0] == "SAVER"):
#             im.saver(
#                 matrix_normalized,
#                 matrix_imputed,
#                 script_path,
#                 rscript_path)
#             print("\tImputed matrix has been written into %s" % matrix_imputed)
#             print("[ %s ] Imputation model done!\n\n" % time.asctime())
#         elif(criteria[0] == "DrImpute"):
#             im.drimpute(
#                 matrix_normalized,
#                 matrix_imputed,
#                 script_path,
#                 rscript_path)
#             print("\tImputed matrix has been written into %s" % matrix_imputed)
#             print("[ %s ] Imputation model done!\n\n" % time.asctime())
#
    return matrix_imputed


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


def cell_trajectory(
        matrix_imputed,
        output_dir,
        parameter_config,
        cluster_file):
    trajectory_parameter = parameter_config['Trajectory']
    rscript_path = parameter_config['Rscript']
    script_path = os.path.join(output_dir, "shell")

    if(trajectory_parameter == 'None'):
        print(
            "[ %s ] Trajectory module does not execute!\n\n" %
            time.asctime())
    else:
        criteria = re.split('\\s', trajectory_parameter)
        print(
            "[ %s ] processing trajectory analysis, the tool and the parameter are as bellow:\n" %
            time.asctime())
        print("\tTools: %s" % criteria[0])
        print("\tParametaers: %s\n" % " ".join(criteria[1:]))
        tr.trajectory_by_fateid(
            matrix_imputed,
            output_dir,
            script_path,
            rscript_path,
            cluster_file)


def go_enrichment(gene_file, output_dir, parameter_config, gene_type=1):
    goa_file = parameter_config['goa']
    obo_file = parameter_config['obo']

    if os.access(goa_file, os.F_OK) and os.access(obo_file, os.F_OK):
        pass
    else:
        print("Warning: goa_file or obo_file cannot be accessed, the default goa and obo files  stored in docs were used instead!")
        current_path = os.path.dirname(__file__)
        goa_file = current_path.strip("core") + "docs/GO/goa_human_isoform.gaf"
        obo_file = current_path.strip("core") + "docs/GO/go.obo"

    rscript_path = parameter_config['Rscript']
    GO_parameter = parameter_config['GO']
    if (GO_parameter == 'None'):
        print("[ %s ] GO module does not execute!\n\n" % time.asctime())
    else:
        if not os.path.exists(output_dir + "/go_enrichment"):
            os.makedirs(output_dir + "/go_enrichment")
        print("[ %s ] Reading %s ..." % (time.asctime(), goa_file))
        output_dir = output_dir + "/go_enrichment"
        er.enrichment(
            gene_file,
            goa_file,
            obo_file,
            output_dir,
            rscript_path,
            gene_type=1)


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


def run(input_file, config_file, output_dir, gene_lst=None):
    initial_module(output_dir, "")
    parameter_config = read_config(config_file)
    matrix_filtered = filter_matrix(input_file, output_dir, parameter_config)
    matrix_normalized = normalize_matrix(
        matrix_filtered, output_dir, parameter_config)

    matrix_imputed = impute_missing_value(
        matrix_normalized, output_dir, parameter_config)

    cluster_file, markers_file, celltype_file = cell_cluster(
        matrix_imputed, output_dir, parameter_config)

    cell_cell_interaction(
        matrix_filtered,
        celltype_file,
        output_dir,
        parameter_config)

    # cell_trajectory(matrix_imputed, output_dir, parameter_config,
    # cluster_file)

    # go_enrichment(markers_file, output_dir, parameter_config, gene_type=1)
