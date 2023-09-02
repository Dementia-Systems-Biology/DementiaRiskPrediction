![licence](https://badgen.net/badge/Licence/MIT/purple)
![status](https://badgen.net/badge/Status/Complete/green)

<h1 align="center">
Blood-based multivariate methylation risk score
   <br>
for cognitive impairment and dementia
</h1>

<p align="center">
<a href="https://github.com/jarnokoetsier/DementiaRiskPrediction/blob/main/README.md#Background">Background</a>
     ·
<a href="https://github.com/jarnokoetsier/DementiaRiskPrediction/blob/main/README.md#Methods">Methods</a>
     ·
<a href="https://github.com/jarnokoetsier/DementiaRiskPrediction/blob/main/README.md#Shiny">Shiny</a>
     ·
<a href="https://github.com/jarnokoetsier/DementiaRiskPrediction/blob/main/README.md#Software">Software</a>
     ·
<a href="https://github.com/jarnokoetsier/DementiaRiskPrediction/blob/main/README.md#Contact">Contact</a>
</p>

<p align="center">
This repository contains the scripts used for the project <i>"Blood-based multivariate methylation risk score for cognitive impairment and dementia"</i>.
</p>
<br>

## Background
As DNA methylation may act as the molecular link between lifestyle/environment and the biological processes governing health and disease, blood-derived DNA methylation data might be utilized for the early identification of persons at risk of developing dementia. Therefore, our research aim is to establish a robust model for predicting a person's midlife dementia risk.

![Methylation](/Images/Methylation.PNG?raw=true "Methylation")


## Methods
The applied methodology consists of four main steps:
1. **Model generation** using the DNA methylation data of the [EXTEND](https://github.com/jarnokoetsier/DementiaRiskPrediction/tree/main/EXTEND) and [EMIF-AD](https://github.com/jarnokoetsier/DementiaRiskPrediction/tree/main/EMIF-AD) cohorts. 
2. **Model validation** in the [EMIF-AD](https://github.com/jarnokoetsier/DementiaRiskPrediction/tree/main/EMIF-AD), [PPMI](https://github.com/jarnokoetsier/DementiaRiskPrediction/tree/main/PPMI), and [ADNI](https://github.com/jarnokoetsier/DementiaRiskPrediction/tree/main/ADNI) cohorts.
3. **Model interpretation** ([here](https://github.com/jarnokoetsier/DementiaRiskPrediction/tree/main/Models/ModelInterpretation)).
4. **Model extension** with genetic and cerebral spinal fluid biomarkers in the [EMIF-AD](https://github.com/jarnokoetsier/DementiaRiskPrediction/tree/main/EMIF-AD) cohort.

## Shiny
You can also generate the multivariate risk scores for cognitive impairment and dementia youself in an R Shiny application:
`runApp('Shiny')`
Please note that it takes a long time to launch the app for the first time, because large models (> 10 GB) are being downloaded from [Zenodo](https://zenodo.org/record/8306113).

## Software
* `R version 4.1.2`
* `RStudio version 2021.09.02+382`

## Contact
Feel free to contact me via email: jarno.koetsier@maastrichtuniversity.nl
