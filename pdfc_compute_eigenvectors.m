function outdata = pdfc_compute_eigenvectors(indata,T,options)
% outdata = pdfc_compute_eigenvectors(indata,T,options)
% 
% Compute leading eigenvector series from indata. This includes computing
% regional Hilbert phase series for every scan session and subsequently
% establish the phase coherence map and its leading eigenvector for every
% time point.
% Input:
%%% indata - all data concatenated into one nxp array
%%% T - Nx1 cell array, where each element represents a subject. Each
% element should be a N_ix1 vector of the length (in samples) of each
% scan session for subject i. E.g., T{1} = [300,300] if the first subject
% has two acquired scan sessions, both of length 300 samples.
%%% options - struct with elements options.P (number of regions) and
% options.flip_eigenvectors if a sign flip should be applied to
% eigenvectors, to ensure the majority of their elements is negative
% Output:
% outdata - data matrix of the same size as indata.
%
% Anders S Olsen April - November 2021, October 2022
% Neurobiology Research Unit, Copenhagen University Hospital Rigshospitalet

T_concat = [T{:}];
sessionstartindices = [1,cumsum(T_concat)+1];
outdata = nan(size(indata));

pb = CmdLineProgressBar('Computing leading eigenvectors for session ');


for ses = 1:length(T_concat)
    
    idx_ses = sessionstartindices(ses):sessionstartindices(ses+1)-1;
    indatases = indata(idx_ses,:);
    
    HilbertMatrix = angle(hilbert(indatases)); % Compute BOLD phases
    
%     disp(['Computing leading eigenvectors for session ',num2str(ses)])
    pb.print(ses,length(T_concat))
    
    for tt = 1:T_concat(ses) %num timepoints in session
        
        % Fill in coherence matrix
        cohmat = cos(HilbertMatrix(tt,:)'-HilbertMatrix(tt,:));
        
        % Get the leading eigenvector
        [V1,~]=eigs(cohmat,1);
        
        if options.flip_eigenvectors
            
            % Make sure the largest component is negative
            if mean(V1>0)>.5
                V1=-V1;
            end
        end
        
        outdata(idx_ses(tt),:) = V1;
        
    end
    
end
disp('Done computing leading eigenvectors')

