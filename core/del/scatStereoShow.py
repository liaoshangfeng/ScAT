import os
import argparse

parser = argparse.ArgumentParser(
    prog='SCAT StereoShow', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('StereoShow', help='Perform cell clustering by Seurat with Stereo-seq data')
parser.add_argument('--input',dest='input',type=str,required=True,
        help='input tsv filename')

parser.add_argument('--binsize', dest='binsize',type=int,required=True,
        help='Input seurate rds spatial transcriptome data')


parser.add_argument('--sample',dest='sample',type=str,required=True,
        help='sample ID, will be used as output prefix and seurat object ident')

parser.add_argument('--tissue',dest='tissue',type=str,default=None,
        help='csv format file listed cell ID which can be used to lasso, should in x_y format')


parser.add_argument('--outdir',dest='outdir',type=str,required=True,
        help='directory where to save the output files, all output files will be indexed by sample ID')

parser.add_argument('--Rscript',dest='rscript',type=str,required=False,
        default='/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/R411/bin/Rscript',help='The path of Rscript')


parser.add_argument('-r', '--rds', dest = 'rds', help = 'save the pre-seurat matrix in rds format')
parser.add_argument('-5', '--h5ad', dest = 'h5ad', help = 'save the pre-seurat matrix in h5ad format')

parser.add_argument('--minCount', dest = 'minCount', default = 0, type = int, help = 'minimum UMI number')
parser.add_argument('--maxCount', dest = 'maxCount', type = int, help = 'maximum UMI number')
parser.add_argument('--minFeature', dest = 'minFeature', default = 0, type = int, help = 'minimum Feature number')
parser.add_argument('--maxFeature', dest = 'maxFeature', type = int, help = 'maximum Feature number')
parser.add_argument('--vg', dest = 'vg', default = 3000, type = int, help = 'number of variable genes, default 3000')
parser.add_argument('--pc', dest = 'pc', default = 30, type = int, help = 'number of PC to use, default 30')
parser.add_argument('--resolution', dest = 'resolution', required=True,default = 0.5, help = 'cluster resolution, default 0.5')

parser.add_argument('--pointSize', dest = 'pointSize', default = 0.2, help = 'point size of spatial plot, default 0.2')
parser.add_argument('--colors', dest = 'colors', default = 70, type = int, help = 'colors palette, one of c(25, 70), default 70')
args = parser.parse_args()


def r_code_filepath():
    # filepath=os.path.join(os.path.split(os.path.realpath(__file__))[0],'R/Stereo_seurat.zr.R')
    filepath=os.path.join(os.path.split(os.path.realpath(__file__))[0],'R/Stereo_seurat.zr_v2.R')
    return filepath

def main():
    if os.path.exists(args.outdir)==False:
        os.system('mkdir '+args.outdir)
    r_script=r_code_filepath()
    if args.maxFeature==None and args.maxCount==None:
        order= args.rscript +' '+r_script +' --input  '+args.input +' --binsize '+str(args.binsize) +' --sample '+args.sample  + \
              ' --out '+args.outdir +' --minCount ' +str(args.minCount) +' --minFeature '+str(args.minFeature)  +' --vg '+str(args.vg)+' --pc '+str(args.pc) +' --resolution '+str(args.resolution) + \
              ' --pointSize '+str(args.pointSize) +' --colors '+str(args.colors)
    else:
       order= args.rscript +' '+r_script +' --input  '+args.input +' --binsize '+str(args.binsize) +' --sample '+args.sample  + \
              ' --out '+args.outdir +' --minCount ' +str(args.minCount) +' --maxCount ' +str(args.maxCount) +' --minFeature '+str(args.minFeature)  +' --maxFeature '+str(args.maxFeature) + ' --vg '+str(args.vg)+' --pc '+str(args.pc) +' --resolution '+str(args.resolution) + \
              ' --pointSize '+str(args.pointSize) +' --colors '+str(args.colors)
 
    print(order)
    os.system(order)


if __name__ == "__main__":
    main()



