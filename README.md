# Psilocybin dFC

This repository contains a set of matlab and R functions to uncover connectivity structures in data using [Leading Eigenvector Dynamics Analysis](https://sites.google.com/site/cvjoanacabral/codes/leida-leading-eigenvector-dynamics-analysis) and diametrical clustering. These functions were used to generate the results in:

>[Psilocybin modulation of dynamic functional connectivity is associated with plasma psilocin and subjective effects ](https://nodesource.com/products/nsolid) 
Anders S. Olsen, Anders Lykkebo-Valløe, Brice Ozenne, Martin K. Madsen, Dea S. Stenbæk, Sophia Armand, Morten Mørup, Melanie Ganz, Gitte M. Knudsen, Patrick M. Fisher.

### Format of inputs
*indata* - all data (post preprocessing) concatenated into one nxp double array, where p is the number of regions (data dimensionality).

*T* - Nx1 cell array, where each element represents a subject. Each element in *T* should be a vector of the length (in samples) of each scan session for subject i. E.g., T{1} = [300,300] if the first subject has two acquired scan sessions, both of length 300 samples.

*covariates* - struct with each covariate being its own field. Each field should contain a Nx1 cell array of the same format as T with the associated covariate level for each scan session. E.g.,  covariate.PPL{1} = [1.93, 2.85] for the first subject.

*options* - struct specifying parameters for the analysis workflow. Please see *pdfc_main* for an example setup. Mandatory fields are: 
*options.min_k* (int), *options.max_k* (int), *options.run_diam* (0/1), and *options.run_kmeans* (0/1). 
Optional fields are: 
*options.flip_eigenvectors*[0], *options.kmeansIterMax*[1000], *options.kmeansRepl*[5], *options.kmeansInit*['plus'], *options.run_frac_occ*[0], *options.run_dwell_time*[0], *options.run_perm_FO*[0], *options.DTintervals*[min,mean,max], *options.seed*[empty], *options.ROIlabels*[cell array 1:p]


### Overview of functions
*pdfc_main.m*
Main script to specify options and data for the subsequent functions. Please see this script to get an overview of data format requirements.

*pdfc_check_input.m*
A collection of failsafes and checks for the input data and options structure.

*pdfc_compute_eigenvectors.m*
Establishes regional Hilbert phase series for every scan session and subsequently the leading eigenvector of the phase coherence map for every timepoint.

*pdfc_cluster_extractmeasures.m*
Clusters the concatenated leading eigenvectors using either diametrical clustering or Euclidean k-means with an eigenvector sign-flip procedure (needs to be specified in *options*). Subsequently, fractional occurrence and dwell time is computed for every estimated cluster centroid. 

*pdfc_diametrical_clustering.m*
K-means procedure implementing algorithm 2 in [The multivariate Watson distribution: Maximum-likelihood estimation and other aspects](https://www.sciencedirect.com/science/article/pii/S0047259X12002084) (Sra S, Karp D, 2012). Cluster centroids may either be initialized using a uniform distribution or using k-means++.

*pdfc_diametrical_clustering_plusplus.m*
A K-means++ procedure modified for diametrical clustering for optimal cluster centroid initialization.

*pdfc_analyzeclusteringdata.m*
Statistical evaluations of fractional occurrence and dwell time. Fractional occurrence is analyzed for every centroid using a linear mixed-effects model. It is possible to run permutation testing with max-T correction to control FWER (see below). Dwell time is analyzed using Cox proportional hazards survival analysis with a frailty element for subject (equivalent to random effect). Results are collected in an output matlab cell with statistical estimates for every centroid. For permutation testing and survival analysis, csv files are saved locally and R functions are called. 

*pdfc_permmaxT.R*
Max-T correction and permutation testing for the fixed effect in a linear mixed-effects model. For permutation testing, *pdfc_permlme_function.R* (below) is called. 

*pdfc_permlme_function.R*
An implementation of [Permutation Tests for Random Effects in Linear Mixed Models](https://www.jstor.org/stable/23270450?seq=1#metadata_info_tab_contents) (Lee OE, Braun TM 2012)

*pdfc_cox_frailty.R*
Cox proportional hazards frailty model using the "survival" package in R. 


## License
Neurobiology Research Unit, Copenhagen University Hospital Rigshospitalet

Free to use but please cite our paper if you do.

## Bugs and improvements
Please direct potential bugs or suggestions for improvement to anders.s.olsen@nru.dk

Several functions for visualization etc. are planned to be uploaded in 2022
