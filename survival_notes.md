# Methods to find best cutoff for survival

- `KMplotter` - the Kaplan Meier plotter is capable to assess the effect of 54,675 genes on survival using 18,674 cancer samples. These include 5,143 breast, 1,816 ovarian, 2,437 lung, 364 liver, 1,065 gastric cancer patients with relapse-free and overall survival data. The miRNA subsystems include additional 11,456 samples from 20 different cancer types. Primary purpose of the tool is a meta-analysis based biomarker assessment.
    - Györffy, Balazs, Andras Lanczky, Aron C. Eklund, Carsten Denkert, Jan Budczies, Qiyuan Li, and Zoltan Szallasi. “An Online Survival Analysis Tool to Rapidly Assess the Effect of 22,277 Genes on Breast Cancer Prognosis Using Microarray Data of 1,809 Patients.” Breast Cancer Research and Treatment 123, no. 3 (October 2010): 725–31. https://doi.org/10.1007/s10549-009-0674-9. - cutoff selection for survival by scanning gene expression range.

- `ctree` function for automatic cutoff finding and building a regression tree out of multiple covariates. `partykit::ctree()`. 
    - Hothorn, Torsten, Kurt Hornik, and Achim Zeileis. “Ctree: Conditional Inference Trees.” The Comprehensive R Archive Network, 2015, 1–34.

- `Cutoff Finder` - web tool for finding optimal dichotomization with respect to an outcome or survival variable. Five methods. http://molpath.charite.de/cutoff/
    - Budczies, Jan, Frederick Klauschen, Bruno V. Sinn, Balázs Győrffy, Wolfgang D. Schmitt, Silvia Darb-Esfahani, and Carsten Denkert. “Cutoff Finder: A Comprehensive and Straightforward Web Application Enabling Rapid Biomarker Cutoff Optimization.” PloS One 7, no. 12 (2012): e51862. https://doi.org/10.1371/journal.pone.0051862.
