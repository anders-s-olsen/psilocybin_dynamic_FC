function [options] = pdfc_check_input(indata,T,covariates,options)

% Check input into the options structure


if sum([T{:}])~=size(indata,1)
    error('T is misspecified (number of samples does not match indata)')
    return
end
% dependencies and failsafe
if options.run_kmeans && options.run_diam
    error('ERROR: Choose only one clustering method. Either, options.run_diam or options.run_kmeans')
    return
end

if ~options.run_diam && ~options.run_kmeans
    error('ERROR: Please choose a clustering method')
    return
end

if options.run_diam && (options.flip_eigenvectors)
    disp('Diametrical clustering not compatible with sign flipping, ignoring sign flip')
    options.flip_eigenvectors = 0;
end

options.P = size(indata,2);

if ~isfield(options,'seed')
    options.seed = [];
end

if ~isempty(covariates)
    options.covnames = fieldnames(covariates);
    options.numcovs = length(options.covnames);
    
    if ~isfield(options,'plot')||(isfield(options,'plot')&&isempty(options.plot))
        options.plot.frac_occ = zeros(1,options.numcovs);
        options.plot.dwell_time = zeros(1,options.numcovs);
        options.plot.entropy = zeros(1,options.numcovs);
        options.plot.centroids = 0;
    else
        if ~isfield(options.plot,'frac_occ')||isempty(options.plot.frac_occ)
            options.plot.frac_occ = zeros(1,options.numcovs);
        else
            if length(options.plot.frac_occ)~=options.numcovs
                options.plot.frac_occ = repelem(options.plot.frac_occ(1),options.numcovs);
            end
        end
        if ~isfield(options.plot,'dwell_time')||isempty(options.plot.dwell_time)
            options.plot.dwell_time = zeros(1,options.numcovs);
        else
            if length(options.plot.dwell_time)~=options.numcovs
                options.plot.dwell_time = repelem(options.plot.dwell_time(1),options.numcovs);
            end
        end
        if ~isfield(options,'run_frac_occ')||isempty(options.run_frac_occ)
            options.run_frac_occ = zeros(1,options.numcovs);
        else
            if length(options.run_frac_occ)~=options.numcovs
                options.run_frac_occ = repelem(options.run_frac_occ(1),options.numcovs);
            end
        end
        if ~isfield(options,'run_dwell_time')||isempty(options.run_dwell_time)
            options.run_dwell_time = zeros(1,options.numcovs);
        else
            if length(options.run_dwell_time)~=options.numcovs
                options.run_dwell_time = repelem(options.run_dwell_time(1),options.numcovs);
            end
        end
        if ~isfield(options,'run_perm_FO')||isempty(options.run_perm_FO)
            options.run_perm_FO = zeros(1,options.numcovs);
        else
            if length(options.run_perm_FO)~=options.numcovs
                options.run_perm_FO = repelem(options.run_perm_FO(1),options.numcovs);
            end
        end
    end
    
    
    for cov = 1:options.numcovs
        if size([covariates.(options.covnames{cov}){:}])~=size([T{:}])
            error(['Wrong size of covariate input: ',options.covnames{cov}])
            return
        end
    end
    
    
else
    error('No covariates specified.')
    
end


if any(options.run_dwell_time) && ...
        (~isfield(options,'DTintervals')||isempty(options.DTintervals))
    for cov = 1:options.numcovs
        options.DTintervals(cov,1) = min([covariates.(options.covnames{cov}){:}]);
        options.DTintervals(cov,2) = mean([covariates.(options.covnames{cov}){:}]);
        options.DTintervals(cov,3) = max([covariates.(options.covnames{cov}){:}]);
    end

elseif any(options.run_dwell_time) && any(options.plot.dwell_time) && ...
        size(options.DTintervals,1)~=options.numcovs
    disp('Please specify covariate increments in integers for all covariates')
end

if ~isfield(options,'alpha')||isempty(options.alpha)
    options.alpha = 0.05;
end

if ~isfield(options,'ROIlabels')||isempty(options.ROIlabels)
    ROIlabels = {};
    for ROI = 1:options.P
        ROIlabels{ROI} = ['ROI_',num2str(ROI)];
    end
end

if ~isfield(options,'kmeansIterMax')||isempty(options.kmeansIterMax)
    options.kmeansIterMax = 500;
end
if ~isfield(options,'kmeansRepl')||isempty(options.kmeansRepl)
    options.kmeansRepl = 5;
end
if ~isfield(options,'min_k')||isempty(options.min_k)
    disp('Minimum number of states not set, computing model for min_k=2')
    options.min_k = 2;
end
if ~isfield(options,'max_k')||isempty(options.max_k)
    disp('Maximum number of states not set, computing model for max_k=2')
    options.max_k = 12;
end
if ~isfield(options,'run_perm_FO')||isempty(options.run_perm_FO)
    options.run_perm_FO = 0;
end
if ~isfield(options,'run_dwell_time')||isempty(options.run_dwell_time)
    options.run_dwell_time = 1;
end
if ~isfield(options,'run_frac_occ')||isempty(options.run_frac_occ)
    options.run_frac_occ = 1;
end
if any(options.plot.frac_occ) && ~any(options.run_frac_occ)
    disp('Cannot plot fractional occurrence without running stats')
end
if any(options.plot.dwell_time) && ~any(options.run_dwell_time)
    disp('Cannot plot dwell time without running stats')
end






