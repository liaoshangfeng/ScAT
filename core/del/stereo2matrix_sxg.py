#!/usr/bin/env python
# coding: utf-8

# This script can transform the stereo matrix to normal matrix (spot x gene)

import os
import sys
import datetime
import gzip


def transform_matrix(input_file, output_dir):
    output_file = os.path.join(output_dir, ("").join(
        [os.path.basename(input_file).strip('.txt'), '_matrix.txt.gz']))

    print("output is %s \n" % output_file)

    spot_mark = {}
    gene_count = {}

    begin_time = datetime.datetime.now()
    print(
        "\n[%s] The transformation is working \n" %
        (datetime.datetime.now()))

    with open(input_file, 'r') as f_in:
        for line in f_in.readlines()[1:]:
            content = line.strip().split("\t")
            if (len(content) < 4):
                print(content)
                continue
            gene_id = content[0]
            spot_id = content[1] + 'x' + content[2]
            count = content[3]
            if spot_id not in spot_mark:
                spot_mark[spot_id] = ''
            if gene_id not in gene_count:
                gene_count[gene_id] = {}
            gene_count[gene_id][spot_id] = count

    print("\n[%s] File has been read \n" % datetime.datetime.now())
    all_spot = [i for i in spot_mark.keys()]
    print("Number of spot is %s.\n" % len(all_spot))
    print("Number of gene is %s.\n" % len(gene_count))

    with gzip.open(output_file, 'wb') as f_out:
        headline = '\t'.join(all_spot) + '\n'
        f_out.write(headline.encode('utf-8'))
        for gene_id in gene_count.keys():

            spot_exp = dict.fromkeys(all_spot, 0)
            spot_count = gene_count[gene_id]
            spot_exp.update(spot_count)
            exp = [str(spot_exp[i]) for i in all_spot]
            exp_line = gene_id + '\t' + '\t'.join(exp) + '\n'
            f_out.write(exp_line.encode('utf-8'))

    print(
        "\n[%s] The matrix has been transformed and stored in file: %s\n" %
        (datetime.datetime.now(), output_file))
    end_time = datetime.datetime.now()
    run_time = end_time - begin_time
    print("Script running time is %s.\n" % run_time)


if __name__ == '__main__':
    transform_matrix(sys.argv[1], sys.argv[2])
