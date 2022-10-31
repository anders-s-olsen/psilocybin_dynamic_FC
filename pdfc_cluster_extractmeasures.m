function [resultstbl,survtbl,C] = pdfc_cluster_extractmeasures(eigenvectors,T,covariates,options)
% [resultstbl,survtbl,C] = pdfc_cluster_extractmeasures(eigenvectors,T,covariates,options)
% Cluster leading eigenvectors according to a range of k-values specified
% in options.min_k and options.max_k. Cluster centroids, their fractional
% occurrences and dwell times are collected in the two outputs resultscell
% and survtbl.
% Input:
%%% eigenvectors - all eigenvectors concatenated into one nxp array
%%% T - Nx1 cell array, where each element represents a subject. Each
% element should be a N_ix1 vector of the length (in samples) of each
% scan session for subject i. E.g., T{1} = [300,300] if the first subject
% has two acquired scan sessions, both of length 300 samples.
%%% covariates - struct with each covariate being its own field. Each field
% should contain a Nx1 cell array of the same format as T with the
% associated covariate level for each scan session. E.g.,
% covariate.PPL{1} = [1.93, 2.85] for the first subject.
%%% options - struct with many fields as specified in pdfc_check_input.m
% Output:
%%% resultscell - Cell array with information on each scan session and
% centroid, including fractional occurrences and covariate levels
%%% survtbl - table with information on all state occurrences including the
% number of active samples in row.
%%% C - Cell array wit computed centroids
%
% Anders S Olsen April - November 2021, October 2022
% Neurobiology Research Unit, Copenhagen University Hospital Rigshospitalet

% initialize results table
resultstbl_variablenames = [{'Subject','Session'},options.covnames',...
    {'N_centroids','Current_centroid','Fractional_occupancy'}];
resultstbl_variableclass = [{'int16','int16'},options.covclass,...
    {'int16','int16','double'}];

resultstbl  = table('Size',[0,numel(resultstbl_variablenames)],...
    'VariableNames', resultstbl_variablenames,...
    'VariableTypes',resultstbl_variableclass);


% Initialize survtbl
survtbl_variablenames = [{'Subject','Session'},options.covnames',...
    {'N_centroids','start_state','end_state','start_time_s','end_time_s','event'}];
survtbl_variableclass = [{'int16','int16'},options.covclass,...
    {'int16','int16','int16','double','double','int16'}];
survtbl = table('Size',[0,numel(survtbl_variablenames)],...
    'VariableNames', survtbl_variablenames,...
    'VariableTypes',survtbl_variableclass);

T_concat = [T{:}];
sessionstartindices = [1,cumsum(T_concat)+1];
tblcount = 1;


for k = options.min_k:options.max_k
    disp(['Clustering: k = ',num2str(k)])
    
    if options.run_diam
        
        [idx, C{k}] = pdfc_diametrical_clustering(eigenvectors,k,options.kmeansIterMax,options.kmeansRepl,options.kmeansInit,options.seed,options.parallel);
        
    elseif options.run_kmeans
        
        [idx, C{k}] = kmeans(eigenvectors,k,'MaxIter',options.kmeansIterMax,'Replicates',options.kmeansRepl,'Start',options.kmeansInit);
        
    end
    
    if any([options.run_frac_occ,options.run_dwell_time])
    disp('Computing summary measures')
    else
        continue
    end
    
    sescounter = 1;
    
    for sub = 1:length(T)
        
        
        for ses = 1:length(T{sub})     % loop through sessions
            
            % extract state sequence for this subject and session
            stateseq = idx(sessionstartindices(sescounter):sessionstartindices(sescounter+1)-1);
            
            % Loop through cluster centroids to compute fractional occ
            if any(options.run_frac_occ)
                
                h = height(resultstbl)+1; %start row
                resultstbl.Subject(h:h+k-1) = sub;
                resultstbl.Session(h:h+k-1) = ses;
                resultstbl.N_centroids(h:h+k-1) = k;
                for elm = 1:options.numcovs
                    if ~isempty(covariates.(options.covnames{elm}){sub})
                        resultstbl.(options.covnames{elm})(h:h+k-1) = covariates.(options.covnames{elm}){sub}(ses);
                    end
                end
                
                for centroid = 1:k
                    % compute fractional occupancy
                    frac_occ = sum(stateseq==centroid)/numel(stateseq);
                    
                    % Fill in output cell for this centroid
                    resultstbl.Current_centroid(h+centroid-1) = centroid;
                    resultstbl.Fractional_occupancy(h+centroid-1) = frac_occ;
                    
                end
            end
            
            if any(options.run_dwell_time)
                % dwell time information. We locate all state transition
                % indices and collect information in separate rows of
                % survtbl, including covariate levels.
                
                h = height(survtbl)+1; %start row
                f_diff = find(diff(stateseq)~=0);
                % We remove first state in session, since we have no
                % information on the true length of preceding
                % activation time
                % perform right censoring for the last state in session

                survtbl.start_time_s(h:h+numel(f_diff)) = options.TR * [0;f_diff];
                survtbl.end_time_s(h:h+numel(f_diff))   = options.TR * [f_diff;length(stateseq)];
                survtbl.start_state(h:h+numel(f_diff))  = [stateseq(f_diff);stateseq(f_diff(end)+1)];
                survtbl.end_state(h:h+numel(f_diff))    = [stateseq(f_diff+1);0];
                survtbl.event(h:h+numel(f_diff))        = [0;ones(numel(f_diff)-1,1);0];
                
                % input standard information into survtbl
                survtbl.Subject(h:h+numel(f_diff))  = sub;
                survtbl.Session(h:h+numel(f_diff))  = ses;
                survtbl.N_centroids(h:h+numel(f_diff)) = k;
                for elm = 1:options.numcovs
                    if ~isempty(covariates.(options.covnames{elm}){sub})
                        survtbl.(options.covnames{elm})(h:h+numel(f_diff)) = covariates.(options.covnames{elm}){sub}(ses);
                    end
                end
                
                
            end
            
            sescounter = sescounter +1;
        end %ses
        
        
    end %sub
end %k
% toc
fprintf('Done with clustering!\n')