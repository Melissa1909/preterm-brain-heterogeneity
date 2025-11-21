# <a name="top"></a>Heterogeneous, temporally consistent, and plastic brain development after preterm birth

This repository accompanies the paper:

Melissa Thalhammer, Jakob Seidlitz, Antonia Neubauer, Aurore Menegaux, Benita Schmitz-Koep, Maria A. Di Biase, Julia Schulz, Lena Dorfschmidt, Richard A. I. Bethlehem, Aaron Alexander-Bloch, Chris Adamson, Gareth Ball, Joana Sa de Almeida, Richard Beare, Claus Zimmer, Marcel Daamen, Henning Boecker, Peter Bartmann, Dieter Wolke, Dennis M. Hedderich, and Christian Sorg (2024). *Heterogeneous, temporally consistent, and plastic brain development after preterm birth.*

[![Nature Communications DOI](https://img.shields.io/static/v1?label=Nature%20Communications&message=doi%3A10.1038%2Fs41467-025-63967-1&color=E63323&style=flat-square)](https://doi.org/10.1038/s41467-025-63967-1)
[![biorxiv](https://img.shields.io/badge/biorxiv-10.1101%2F2024.12.06.627134-blue?style=flat)](https://doi.org/10.1101/2024.12.06.627134)


## Abstract
The current view of neurodevelopment after preterm birth presents a strong paradox: diverse neurocognitive outcomes suggest heterogeneous neurodevelopment, yet numerous brain imaging studies focusing on average dysmaturation imply largely uniform aberrations across individuals. Here we show both, spatially heterogeneous individual brain abnormality patterns (IBAPs) but with consistent underlying biological mechanisms of injury and plasticity. Using cross-sectional structural magnetic resonance imaging data from preterm neonates and longitudinal data from preterm children and adults in a normative reference framework, we demonstrate that brain development after preterm birth is highly heterogeneous in both severity and patterns of deviations. Individual brain abnormalities were also consistent for their extent and location along the life course, associated with glial cell underpinnings, and plastic for influences of the early social environment. Thus, IBAPs of preterm birth are spatially heterogenous, temporally consistent for extent, spatial location, and cellular underpinnings, and plastic for social-environmental impacts. Our findings extend conventional views of preterm neurodevelopment, revealing a nuanced landscape of individual variation, with consistent commonalities between subjects. This integrated perspective of preterm neurodevelopment implies more targeted theranostic intervention strategies, specifically integrating brain charts and imaging at birth, as well as social interventions during early development.


## License and citation
[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey)](http://creativecommons.org/licenses/by-nc-sa/4.0/)  

Content of this repository is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-nc-sa/4.0/).

Please cite the publication when using the data or code provided in this repository. Third-party data included in this repository requires direct citation if used in other contexts. This applies to all MRI imaging data, normative models, or gene expression data. See Methods and Supplementary Information of the publication for details and resources.


## Installation
This repository contains the code to reproduce the analyses and figures of the publication. The code is written mainly in Python and R, only Fig. 6 is generated using MATLAB. The installation steps are described below. Analyses were carried out on a MacBook Air with macOS Sonoma 14.6.1 and might require adjustments for other operating systems, especially Windows. The code was additionally tested on a HP Pro notebook with Windows 11 Pro. The installation command differ slightly for Windows, please individually adjust the commands accordingly.

### Software requirements
- [Python 3.11.9](https://www.python.org/downloads/release/python-3119/)
- [R 4.4.1](https://cran.r-project.org/)
- [MATLAB 2023b](https://de.mathworks.com/help/install/ug/install-products-with-internet-connection.html) (for Fig. 10 only)

### Installation steps
*Expected time: 10 minutes*
1. Create an empty directory and navigate into it:
    ```bash
    mkdir /path/to/your/setup && cd /path/to/your/setup
    ```
2. Clone this repository (including submodules `Lifespan` containing the pretrained normative models and code to refit them) into the directory:
    ```bash
    git clone --recursive https://github.com/Melissa1909/preterm-brain-heterogeneity 
    ```
3. Create a virtual environment with Python 3.11.9 and activate it:
    ```bash 
    python3.11 -m venv .venv 
    source .venv/bin/activate
    ```
    On Windows, activation is done with:
        ```bash
        .venv\Scripts\activate
        ```
4. Install python dependencies:
    ```bash
    pip install -r requirements.txt
    ```
5. Install the required R packages (run this in an R environment):
    ```R
    install.packages(c("ggplot2", "gamlss", "tidyverse", "devtools", "invgamma"))
    devtools::install_github("jcbeer/longCombat")
    ```
    Furthermore, install the PROCESS toolbox for R. You can download it [here](https://www.afhayes.com/public/processv43.zip), unzip it, and place it into `/code` as well.
6. Install the required MATLAB package (for Fig. 6 only):
    ```bash
    cd preterm-brain-heterogeneity/code
    git clone https://github.com/dutchconnectomelab/Simple-Brain-Plot 
    ```


## Usage
### `Lifespan` submodule
The submodule `Lifespan` contains the [pretrained normative models](https://doi.org/10.1038/s41586-022-04554-y) and code to refit them. Credit goes to R. Bethlehem (@rb643) and J. Seidlitz (@jms290) for providing the code and models. The submodule is included in the repository to ensure reproducibility of the results. Furthermore, a minimal adaptation of the code was necessary to fit the models to the data of the publication (i.e., the only adaptation was to not remove all variables before refitting the models, so that the functions can be run with passing variables). 

### `/data`
The subfolders would contain the data used in the publication but are empty due to data protection regulations. 
- `dHCP`: Neonatal imaging data from the developmental Human Connectome Project. The dataset requires access application [here](https://nda.nih.gov/edit_collection.html?id=3955). Furthermore, data need to be parcellated with the [M-CRIB-S parcellation](https://github.com/DevelopmentalImagingMCRI/MCRIBS), which corresponds to the adult Desikan-Killiany atlas.
- `ABCD`: Children imaging data from the Adolescence Brain Cognitive Development Study. Apply to access [here](https://nda.nih.gov/abcd/requestaccess) and download the data from [release 4.0](https://nda.nih.gov/study.html?id=1299). 
- `BLS`: Adult imaging, behavioral, and perinatal data collected in-house. Access can be granted only upon reasonable request from the authors. 
- `cell_types`: Copy of compiled cell-specific gene set list from all available large-scale single-cell studies of the adult human cortex obtained from the raw [Seidlitz et al.](https://doi.org/10.1038/s41467-020-17051-5) [dataset](https://staticcontent.springer.com/esm/art%3A10.1038%2Fs41467-020-17051-5/MediaObjects/41467_2020_17051_MOESM8_ESM.xlsx).
- `synthetic_BLS-26`: Synthetic dataset for demonstration purposes (see below).


### `Jupyter notebooks`
The Jupyter notebooks contain all the analyses code in sequential order to reproduce the results of the publication. The notebooks are organized as follows:
1. `0.1_dataPrep_*.ipynb`: Preparation of meta- and imaging data for the respective dataset, including quality control, averaging across hemispheres, data harmonization (if applicable), preparation of the input file for subsequent normative modeling. *Expected runtime: < 5 minutes per dataset.*
2. `1_spatialHeterogeneity.ipynb`: Contrasting average dysmaturation outcomes (based on group comparisons with two-sample t-tests) and individual brain abnormality patterns (IBAPs) obtained via normative modeling. We show that brain aberrations after preterm birth are spatially heterogeneous between individuals in all three datasets. This script should be run for all datasets and the modalities cortical thickness `CT` and surface area `SA`. *Expected runtime (depending on dataset size and CPU usage): up to 1 hour per dataset and modality*.
3. `2_macrosaleConsistency.ipynb`: Investigating the temporal consistency of individual brain abnormality patterns (IBAPs) across the life course. We demonstrate that the extent and location of brain abnormalities are consistent in individual development. *Expected runtime: < 5 minutes.*
4. `3_microscaleConsistency.ipynb`: Investigating the cellular underpinnings of IBAPs and their consistency across the life course. The analysis shows that glial cells underpin adult IBAPs, which is complementary to postmortem findings of glial damage in preterm neonates. *Expected runtime: up to 2 hours.*
5. `4_developmentalPlasticity_IQ.ipynb`: Assessing the plasticity of IBAPs in response to the early social environment and effects of adult IBAPs on cognitive outcome variability. We show that the early social environment influences the severity of brain abnormalities in preterm individuals and that the severity of IBAPs is linked to IQ in adulthood. *Expected runtime: < 5 minutes.*
6. `Supp_Rutherford.ipynb`: Reproducing the spatial heterogeneity results in the BLS dataset using a different pre-trained normative model by [Rutherford et al. (2022)](https://doi.org/10.7554/eLife.72904). *Expected runtime: about 5-10 minutes.*

### `MATLAB script`
Parts of Figure 6 were generated with `Fig6_concept.m` in MATLAB 2023b using Simple-Brain-Plot. *Expected runtime: < 5 minutes.*

### `/code`
Contains additional functions for analyses, which are imported in the Jupyter notebooks.

### `/outputs`
Contains the results of the analyses, including figures.


## Minimal working example
Due to data sharing restrictions, we cannot provide the original data used in the publication. However, we provide a synthetic dataset in `data/synthetic_BLS-26` to demonstrate the heterogeneity analysis of preterm adults for CT and SA. The synthetic dataset contains the same variables as the original `BLS-26` dataset, but the data were artificially created by random sampling from a Gaussian distribution with the same mean and standard deviation as the original data. The synthetic data were not used in the publication and do not correspond to the results of the publication.

To run the analysis on the synthetic data, run `1_spatialHeterogeneity.ipynb` with `dataset` set to `synthetic_BLS-26`. The results will be saved in the `outputs/spatial_heterogeneity/synthetic_BLS-26` folder. The synthetic data is only used for demonstration purposes and do not correspond to the results of the publication.
For the other analyses, the original data is required.


## Troubleshooting
### Installation issues
- If you encounter issues with installing the [ENIGMA toolbox](https://enigma-toolbox.readthedocs.io/en/latest/), please try to install the toolbox again manually as described in the [ENIGMA documentation](https://enigma-toolbox.readthedocs.io/en/latest/pages/01.install/index.html). Important: `vtk==9.2.6` is required, I've encountered problems with other versions of vtk.
- Make sure Python and R installations are set in the PATH environment variable. If you encounter issues with the R packages, try to install them manually in RStudio or R console.


## Support and contact 
If you run into issues, please open an issue on this repository or [contact me](mailto:melissa.thalhammer@tum.de), including your operating system and a detailed description of the problem. Feedback is appreciated!
