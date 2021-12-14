% Collection of data analysis scripts pertaining:
% "Psilocybin modulation of dynamic functional connectivity is associated 
% with plasma psilocin and subjective effects" by Olsen AS et al., 2021
% 
% 1) Collect and concatenate data
% 2) Compute leading eigenvector series
% 3) Run k-means or diametrical clustering and compute metrics fractional 
% occurrence and/or dwell time
% 4) Run statistics on fractional occurrence and/or dwell time
%
% Scripts for generating visualizations are planned to be added to the
% github repository github.com/anders-s-olsen/psilocybin_dynamic_FC
%
% Anders S Olsen April - November 2021
% Neurobiology Research Unit, Copenhagen University Hospital Rigshospitalet

clear,close all
warning('off','all')
options.functionpath = [pwd,'/']; % full path to pdfc functions

%% Choose which parts of the pipeline to run

% Define parameters
options.kmeansIterMax = 500;
options.kmeansRepl = 5;
options.min_k = 2;
options.max_k = 10;

% Clustering type
options.run_diam = 1;  % Run diametrical clustering (applicable to LEiDA only)
options.run_kmeans = 0;% Run k-means clustering (applicable to LEiDA only)

% Flip eigenvectors? Not applicable to diametrical clustering but should be
% applied when using Euclidean kmeans
options.flip_eigenvectors = 0;

% Statistics - these are logical vectors where each element is associated
% with the respective covariate in the "covariates" variable
options.run_frac_occ = [0,0];
options.run_dwell_time = [0,0];
options.run_perm_FO = [0,0];       % Do permutation testing in R (only FO)

% Dwell-time intervals for which to compute marginal survival curves from
% the Cox frailty PH model - one row for each covariate. Integer only
options.DTintervals = [0,10,20;0,5,10];

options.seed = 0; %leave empty if no seed

%% NRU specific - remember to delete

sourcepath = '/data1/anderssolsen/np2/';
addpath(genpath(sourcepath))
load([sourcepath,'np2_data_conn_SDI']);

c2_hdr = c2(1,:);
c2 = c2(2:end,:); % remove header information

% exclude data
exclude_art_samples = 0; % Exclude volumes flagged by ART
exclude_p377 = 1;   % Exclude two sets from p377 that are pure motion
options.TR = 2;

% Locate data in the data cell (some rows are empty)
idx_data = find(~cellfun('isempty', c2(:,strcmp(c2_hdr,'data'))));

if exclude_p377
    excl = find(strcmp('yes',c2(:,strcmp(c2_hdr,'exclude')))); %sessions excluded
    [~,memberLOC] = ismember(excl,idx_data);
    idx_data(memberLOC)=[];
end

%% concatenate data and make T and covariate cells

% Initialize
indata=[];sescount=1;sub=0;

for ses = 1:length(idx_data)
    clearvars art_flag inp
    
    if sescount == 1
        sub = sub+1;
    end
    
    
    inp = c2{idx_data(ses),strcmp(c2_hdr,'data')}; %raw data
    
    if exclude_art_samples
        art_flag = c2{idx_data(ses),strcmp(c2_hdr,'art_outlier')};
        
        if ~isempty(art_flag)
            inp(:,art_flag) = [];
        end
    end
    
    indata = [indata;inp'];
    
    T{sub}(sescount) = size(inp,2);
    covariates.PPL{sub}(sescount) = 1020*str2double(c2{idx_data(ses),strcmp(c2_hdr,'psilo_conc')});
    covariates.SDI{sub}(sescount) = c2{idx_data(ses),strcmp(c2_hdr,'SDI')};
    
    % If current session is the last one for the subject, restart counter
    sescount = sescount + 1;
    if idx_data(ses)==max(idx_data)
        sescount = 1;
    elseif ~strcmp(c2{idx_data(ses),strcmp(c2_hdr,'subject')},c2{idx_data(ses+1),strcmp(c2_hdr,'subject')})
        sescount = 1;
    end
    
    
end

%% Format of inputs

%%% indata - all data concatenated into one nxp array

%%% T - Nx1 cell array, where each element represents a subject. Each
% element should be a N_ix1 vector of the length (in samples) of each
% scan session for subject i. E.g., T{1} = [300,300] if the first subject
% has two acquired scan sessions, both of length 300 samples.

%%% covariates - struct with each covariate being its own field. Each field
% should contain a Nx1 cell array of the same format as T with the
% associated covariate level for each scan session. E.g., 
% covariate.PPL{1} = [1.93, 2.85] for the first subject.

%% Run algorithms

options = pdfc_check_input(indata,T,covariates,options);

eigenvectors = pdfc_compute_eigenvectors(indata,T,options);

[resultscell,survtbl] = pdfc_cluster_extractmeasures(eigenvectors,T,covariates,options);

p_cell = pdfc_analyzeclusteringdata(resultscell,survtbl,options);


