# multitemporal_sdm_assessment

## GitHub repository for the paper *Assessing multitemporal calibration for species distribution models*

Here we provide the source code for the computational experiment. A minimal example parametrization is provided (see the comments along the code in the main script, `main_script.R`). For the minimal example runing, please dowload the required climate data and auxiliary shapefiles published at the Mendeley Data repository (this [link](https://data.mendeley.com/datasets/ykg27hk766/1)). Locate the climate data at a folder named `clikmate_data` and the auxiliary shapefiles at a `utils` folder, all in the main project folder. Such data files are all you need to run all the steps of our experiment pipeline. 

For the fully reproduction of our results, please dowload the full climate data, provided at DRYAD (this [link](https://datadryad.org/stash/dataset/doi:10.5061/dryad.8kc1v)). After that, decompress all the files and locate it at `climate_data` folder. For your convenience, we also provide an accessory script (see `decompressor.bat` file) for recursive decompression of the climate data files.

The conceptual description of the experiment pipeline is depicted in the following figure.

<!-- ![Experiment pipeline](/assets/fig01_artificial_sps_large.jpg "Experiment pipeline") -->

<div align="center">
<img src=./assets/fig01_artificial_sps_large.jpg width="900">
</div>

Have a happy coding. :abacus:
