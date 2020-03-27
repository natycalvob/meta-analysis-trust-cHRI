# A Meta Analysis on Children’s Trust in Social Robots. 
(Rebecca Stower &amp; Natalia Calvo-Barajas &amp; Ginevra Castellano &amp; Arvid Kappas)

This repository contains the a set of R scripts to run a Multivariate Mixed-Effect Meta-Analytic Regression Model for Meta Analysis on Children's Trust in Social Robots, but also the information about the studies and methods used for data-extraction. 

| File Name | Description |
| --- | --- |
| main_Meta-Analysis-Trust-cHRI.R | Main script containing the four mixed-effects meta-analytic models |
| MRE_model.R | Multivarite model, univariate model, moderator analysis |
| correlation.R | Imputed correlations |
| effectsizes_variance.R | Computing effect sizes estimates (Cohen's d) based on means, SDs, t-value, and F-values |
| heterogeneity.R | Heterogeneity analysis |
| visualizations.R | Meta-Analytic visualizations |
| MetaAnalysis_Trust_cHRI | Collection of studies related to this research, methods, selection criteria, coding book for data extraction.|

*Note 1: We follow the methodology proposed by Standford Metalab to calculate the effect sizes estimates based on the experimental design, variance, and data extration. Please check the MetaLab repository for more details.* </br>
[Github Metalab](https://github.com/langcog/metalab2)

## Requirements

- R <br />
In order to succesfully run the code, you should install the following libraries in your R environment <br/>
- tidyverse <br />
- RCurl <br />
- metafor <br />
- pwr <br />
- dplyr <br />
- Hmisc <br />
- clubSandwich <br />

## Data Set

You can find the dataset for the current Meta Analysis [here.](https://docs.google.com/spreadsheets/d/e/2PACX-1vQeQHaYFl1a8Pm5oz-k2oYyb6IUpJ7NLeSgSo44wWSCsfYbexgxa7i7ZHha5s8wG3jCNr_dwcsoEFut/pub?output=csv)

## Installation 
Do not forget to install the necessary libraries described in the Requirements paragraph.
Clone the repository in the source folder of your workspace.

```
git clone https://github.com/natycalvob/meta-analysis-trust-cHRI.git
```


## Author
* [Natalia Calvo-Barajas](https://github.com/natycalvob), e-mail: natalia.calvo@it.uu.se
* Rebecca Stower, e-mail: r.stower@jacobs-university.de
