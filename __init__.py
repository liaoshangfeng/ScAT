import sys
# import os
import argparse
description = "Description:\n\n" + \
              "SCAT allows you to perform scRNA-seq analysis, \n" \
              "Call ASE from WASP filter bam \n" \
              "aseQTL to find cis-candidate regulatory SNP\n" \
              "For further information, see the help of each subcommand."

parser = argparse.ArgumentParser(description=description, prog='ASEkit')
# version
parser.add_argument(
    '-v',
    '--version',
    action='version',
    version='%(prog)s 1.0.1',
    help="ssssssssssss"
)
# version = '1.0.0'

# args = parser.parse_args()


def help_info():
    help_info = 'SCAT : A python pipline apply for analyzing scRNA-seq data,\n\n' + \
        'Function 1: Oneclick' + '\n' + 'Perform single cell analysis by oneclick,\n\n' \
        'Function 2: Normalization' + '\n' + 'Perform Normalization of raw data,\n\n' \
        'Function 3: Clustering' + '\n' + 'Perform seurat workflow,\n\n' \
        'Function 4: Annotation' + '\n' + 'Perform cell annotation by Singler, \n\n' \
        'Function 5: Trajectory' + '\n' + 'Perform cell trajectory by Monocle3, \n\n' \
        'Function 6: Enrichment' + '\n' + 'Perform Go Enrichment analysis by ClusterProfiler\n\n'
    return (help_info)


def main():
    # Base_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # sys.path.append(Base_DIR)

    if len(sys.argv) == 1:
        print(help_info())
        sys.exit(1)

    if sys.argv[1] == 'Oneclick':
        from core import scatOneclick
        scatOneclick.main()

    if sys.argv[1] == 'Normalization':
        import scatNormalization
        scatNormalization.main()

    if sys.argv[1] == 'Clustering':
        import scatClustering
        scatClustering.main()

    if sys.argv[1] == 'Enrichment':
        import scatEnrichment
        scatEnrichment.main()

    if sys.argv[1] == 'Trajectory':
        import scatTrajectory
        scatTrajectory.main()

    if sys.argv[1] == 'Annotation':
        import scatAnnotation
        scatAnnotation.main()

    # try:
    #    # command for testing
    #    if sys.argv[1] == 'Clustering':
    #        from . import scatClustering
    #        scatClustering.main()
    #    # elif sys.argv[1] == 'test':
    #    #     test_filepath=os.path.join(os.path.split(os.path.realpath(__file__))[0],'example.data')
    #    #     print('test data filepath is: '+test_filepath)
    #    elif sys.argv[1] in ['version', '--version', '-V', '-v']:
    #        print('[SCAT] version: %s' % version)
    #    # command for help information
    #    elif sys.argv[1] in ['help', '--help', '-h']:
    #        print(help_info())
    #    else:
    #        print('Can not know this parameter, please see README of the package')
    # except Exception:
    #     print('Can not know this parameter,please see  README of the package')


if __name__ == '__main__':
    main()
