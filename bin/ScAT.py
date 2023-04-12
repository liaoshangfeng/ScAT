#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Shangfeng Liao
# @E-mail: liaoshangfeng@genomics.cn
# @Date:   2022-12-21 10:56:00
# @Last Modified by:   Shangfeng Liao
# @Last Modified time: 2023-02-23 10:00:00

import os
import sys
import logging
import argparse
from version import __version__

from scatQcControl import QcControlParser, QcControl 
from scatClustering import ClusteringParser, Clustering
from scatSpatialCluster import SpatialClusterParser, SpatialCluster
from scatAnnotation import AnnotationParser, Annotation
from scatAnnotscCATCH import AnnotscCATCHParser, AnnotscCATCH
from scatCellCellInteraction import CellCellInteractionParser, CellCellInteraction
from scatDeconvSpot import DeconvSpotParser, DeconvSpot
from scatTFregulon import TFregulonParser, TFregulon
from scatSpatialPatternGene import SpatialPatternGeneParser, SpatialPatternGene
from scatEnrichment import EnrichmentParser, Enrichment
from scatTrajectory import TrajectoryParser, Trajectory


def main():
    """
    Add main function argument parsers.
    """
    parser = argparse.ArgumentParser(prog='ScAT',description='Stereo-seq and scRNA-seq Analysis Toolkit (ScAT), a user-friendly Python package that provides multiple modules for scRNA-seq and Stereo-seq data analysis.',
                                  epilog='Use %(prog)s {command} -h to get help on individual commands')
    
    # parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=__version__), help = "Print version info.")   
    parser.add_argument("-v", "--version", action = "store_true", help = "print version info.")
    
    subparsers = parser.add_subparsers(dest = "subcommand")
    QcControlParser(subparsers)
    ClusteringParser(subparsers)
    TrajectoryParser(subparsers)
    AnnotationParser(subparsers)
    AnnotscCATCHParser(subparsers)
    EnrichmentParser(subparsers)
    CellCellInteractionParser(subparsers)
    TFregulonParser(subparsers)
    SpatialClusterParser(subparsers)
    DeconvSpotParser(subparsers)
    SpatialPatternGeneParser(subparsers)
    
    
    logging.basicConfig(format="%(levelname)s: %(message)s", stream=sys.stderr)
    args = parser.parse_args()
    
    version = __version__
    
    
    if args.version:
       print(version)
       exit(0)
    elif args.subcommand == "QcControl":
        QcControl(args)
    elif args.subcommand == "Clustering":
        Clustering(args)
    elif args.subcommand == "Trajectory":
        Trajectory(args)
    elif args.subcommand == "Annotation":
        Annotation(args)
    elif args.subcommand == "AnnotscCATCH":
        AnnotscCATCH(args)
    elif args.subcommand == "Enrichment":
        Enrichment(args)
    elif args.subcommand == "CellCellInteraction":
        CellCellInteraction(args)
    elif args.subcommand == "TFregulon":
        TFregulon(args)
    elif args.subcommand == "SpatialCluster":
        SpatialCluster(args)
    elif args.subcommand == "DeconvSpot":
        DeconvSpot(args)
    elif args.subcommand == "SpatialPatternGene":
        SpatialPatternGene(args)
    else:
       parser.print_help()
       exit(1)
    exit(0)

    
if __name__ == "__main__":
    main()      
