# ScAT

[![python >3.8.8](https://img.shields.io/badge/python-3.8.8-brightgreen)](https://www.python.org/) 

### ScAT: A package designed for joint analysis of Single-cell RNA sequencing and Stereo-seq data
Here, we present Stereo-seq and scRNA-seq Analysis Toolkit (ScAT), a user-friendly Python package that provides multiple modules for scRNA-seq and Stereo-seq data analysis. We use a friendly parameter passing method so that skilled users can modify the pipeline by customizing the parameters. Except for the spatial clustering analysis, ScAT also supports pipelines for integrated analysis of scRNA-seq and Stereo-seq data, such as cell type deconvolution. The data preprocessing processes, including QC control, normalization, and clustering of scRNA-seq and Stereo-seq data, were supported by the wrapper functions for ‘Seurat’ R package.


<p align="center" width="100%">
    <img width="50%" src="https://github.com/liaoshangfeng/ScAT/blob/main/example/ScAT_overview.jpg" alt="" title="ScAT Overview"> 
</p>

# Dependences

[![numpy-1.21.3](https://img.shields.io/badge/numpy-1.21.3-red)](https://github.com/numpy/numpy)
[![pandas-1.2.4](https://img.shields.io/badge/pandas-1.2.4-lightgrey)](https://github.com/pandas-dev/pandas)
[![scanpy-1.8.1](https://img.shields.io/badge/scanpy-1.8.1-blue)](https://github.com/theislab/scanpy)


# Usage

```
wget https://github.com/liaoshangfeng/ScAT
cd ScAT/bin
python3 ScAT -h
```

# Environment

```
conda create -n ScAT
wget https://github.com/liaoshangfeng/ScAT/blob/main/env/install_env.sh
bash install_env.sh
```

# Example

1. [Vignette](https://github.com/liaoshangfeng/ScAT/blob/main/example/ScAT_Vignette.ipynb)

2. [Mouse kidney stereo-seq data analysis tutorial](https://github.com/liaoshangfeng/ScAT/blob/main/tutorial/ScAT_mouse_kidney_tutorial_F2.ipynb)

3. [Mouse kidney snRAN-seq data analysis tutorial](https://github.com/liaoshangfeng/ScAT/blob/main/tutorial/ScA_mouse_kidney_sc_tutorial.ipynb)

# Disclaimer

This tool is for research purpose and not approved for clinical use.


# Coypright

This tool is developed in BGI research.

The copyright holder for this project is BGI research.

All rights reserved.
