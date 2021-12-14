function [resultscell,survtbl] = pdfc_cluster_extractmeasures(eigenvectors,T,covariates,options)
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



% Initialize resultscell
resultscell_hdr{1,1} = 'Subject';
resultscell_hdr{1,2} = 'Session';
for elm = 1:options.numcovs
    resultscell_hdr{1,2+elm} = options.covnames{elm};
end
resultscell_hdr{1,2+options.numcovs+1} = 'N_centroids';
resultscell_hdr{1,2+options.numcovs+2} = 'Current_centroid';
resultscell_hdr{1,2+options.numcovs+3} = 'Fractional_occupancy';
resultscell_hdr{1,2+options.numcovs+4} = 'centroid';
resultscell = resultscell_hdr; %keep header information
cellcount = 2;

% Initialize survtbl
survtbl_hdr = [resultscell_hdr(1,1:end-3),'start_state','end_state','start_time','end_time','event'];
survtbl = table('Size',size(survtbl_hdr),'VariableTypes',...
    [{'int16','int16'},repelem({'double'},options.numcovs),{'int16','int16','int16','double','double','int16'}],...
    'VariableNames',survtbl_hdr);


T_concat = [T{:}];
sessionstartindices = [1,cumsum(T_concat)+1];


for k = options.min_k:options.max_k
    disp(['Clustering: k = ',num2str(k)])
    
    if options.run_diam
        
        [idx, C] = pdfc_diametricalkmeans(eigenvectors,k,options.kmeansIterMax,options.kmeansRepl,'++',options);
        
    elseif options.run_kmeans
        
        [idx, C] = kmeans(eigenvectors,k,'MaxIter',options.kmeansIterMax,'Replicates',options.kmeansRepl);
        
    end
    
    disp('Computing summary measures')
    
    sescounter = 1;
    
    for sub = 1:length(T)
        
        
        for ses = 1:length(T{sub})     % loop through sessions
            
            idx_ses = sessionstartindices(sescounter):sessionstartindices(sescounter+1)-1;
            
            % extract state sequence for this subject and session
            stateseq = idx(idx_ses);
            
            % Loop through cluster centroids to compute fractional occ
            if any(options.run_frac_occ)
                for centroid = 1:k
                    
                    % Fill in output cell for this centroid
                    resultscell{cellcount,strcmp(resultscell_hdr,'Subject')}          = num2str(sub);
                    resultscell{cellcount,strcmp(resultscell_hdr,'Session')}          = num2str(ses);
                    for elm = 1:options.numcovs
                        if ~isempty(covariates.(options.covnames{elm}){sub})
                            resultscell{cellcount,strcmp(resultscell_hdr,options.covnames{elm})} = covariates.(options.covnames{elm}){sub}(ses);
                        end
                    end
                    resultscell{cellcount,strcmp(resultscell_hdr,'N_centroids')}      = num2str(k);
                    resultscell{cellcount,strcmp(resultscell_hdr,'Current_centroid')} = num2str(centroid);
                    
                    % compute fractional occupancy
                    frac_occ = sum(stateseq==centroid)/numel(stateseq);
                    resultscell{cellcount,strcmp(resultscell_hdr,'Fractional_occupancy')} = frac_occ;
                    
                    
                    % only once for every ncentroids, we input cluster
                    % centroids in our results cell
                    if sub==1&&ses==1
                        resultscell{cellcount,strcmp(resultscell_hdr,'centroid')} = C(centroid,:);
                    end
                    cellcount = cellcount + 1;
                    
                end
            end
            
            if any(options.run_dwell_time)
                % dwell time information. We locate all state transition
                % indices and collect information in separate rows of
                % survtbl, including covariate levels.
                
                f_diff = find(diff(stateseq)~=0);
                for trans = 1:numel(f_diff)
                    % We remove first state in session, since we have no
                    % information on the true length of preceding
                    % activation time
                    if trans == 1
                        continue
                    end
                    
                    % input information into survtbl
                    l = height(survtbl);
                    survtbl.Subject(l+1)  = sub;
                    survtbl.Session(l+1)  = ses;
                    for elm = 1:options.numcovs
                        if ~isempty(covariates.(options.covnames{elm}){sub})
                            survtbl.(options.covnames{elm})(l+1) = covariates.(options.covnames{elm}){sub}(ses);
                        end
                    end
                    survtbl.N_centroids(l+1) = k;
                    
                    survtbl.start_time(l+1) = options.TR * f_diff(trans-1);
                    survtbl.end_time(l+1) = options.TR * f_diff(trans);
                    survtbl.start_state(l+1) = stateseq(f_diff(trans));
                    survtbl.end_state(l+1)   = stateseq(f_diff(trans)+1);
                    survtbl.event(l+1)       = 1; % event = state death
                    
                end
                % perform right censoring for the last state in session
                survtbl.start_time(l+1)  = options.TR * f_diff(end);
                survtbl.end_time(l+1)  = options.TR * length(stateseq);
                survtbl.start_state(l+1) = stateseq(f_diff(end)+1);
                survtbl.end_state(l+1)   = NaN;
                survtbl.event(l+1)       = 0; % event = right censored
                
            end 
            
            sescounter = sescounter +1;
        end %ses
        
        
    end %sub
end %k

survtbl(1,:) = [];
% toc
fprintf('Done with clustering!\n')
