% Collection of data analysis scripts pertaining:
% "Psilocybin modulation of dynamic functional connectivity is associated 
% with plasma psilocin and subjective effects" by Olsen AS et al., 2021
% 
% 1) Collect and concatenate data
% 2) Compute leading eigenvector series
% 3) Run k-means or diametrical clustering and compute metrics fractional 
% occurrence and/or dwell time
% 4) Run statistics on fractional occurrence and/or dwell time associations
% with covariates, e.g., plasma psilocin level
%
% Anders S Olsen April - November 2021, October 2022
% Neurobiology Research Unit, Copenhagen University Hospital Rigshospitalet

clear,close all
warning('off','all')
options.functionpath = [pwd,'/']; % full path to pdfc functions

%% Choose which parts of the pipeline to run

% Define parameters
options.kmeansIterMax = 500;
options.kmeansRepl = 5;
options.kmeansInit = 'plus';
options.min_k = 2;
options.max_k = 2;
options.TR = 2;

% Clustering type
options.run_diam = 1;  % Run diametrical clustering
options.run_kmeans = 0;% Run k-means clustering

% Flip eigenvectors? Not applicable to diametrical clustering but should be
% applied when using Euclidean kmeans
options.flip_eigenvectors = 0;

% Statistics - these are logical vectors where each element is associated
% with the respective covariate in the "covariates" variable
options.run_frac_occ = [true,true];
options.run_dwell_time = [true,true];
options.run_perm_FO = [true,true];% Do permutation testing in R (only FO)

% Dwell-time intervals for which to compute marginal survival curves from
% the Cox frailty PH model - one row for each covariate. Integer only
options.DTintervals = [0,10,20;0,5,10];

options.seed = []; %leave empty if no seed

%% Format of inputs

%%% indata - all data concatenated into one nxp array

%%% T - Nx1 cell array, where each element represents a subject. Each
% element should be a N_ix1 vector of the length (in samples) of each
% scan session for subject i. E.g., T{1} = [300,300] if the first subject
% has two acquired scan sessions, both of length 300 samples.

%%% covariates - struct with each covariate being its own field. Each field
% should contain a Nx1 cell array of the same format as T with the
% associated covariate level for each scan session. E.g., 
% covariates.PPL{1} = [1.93, 2.85] for the first subject.

% example using random data:
nsubs = 10;
indata = randn(nsubs*300,90);
T = repelem({300},nsubs);
covariates.PPL = {1.93,2.85,0,1.77,0,0.45,0,2.77,0,3.89};

%% Run algorithms

options = pdfc_check_input(indata,T,covariates,options);

eigenvectors = pdfc_compute_eigenvectors(indata,T,options);

[resultstbl,survtbl] = pdfc_cluster_extractmeasures(eigenvectors,T,covariates,options);

p_cell = pdfc_analyzeclusteringdata(resultstbl,survtbl,options);


