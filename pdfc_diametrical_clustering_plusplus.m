function C = pdfc_diametrical_clustering_plusplus(X,k,stream)

[n,~] = size(X);

% 1) Choose first centroid at random between the samples
% next = randi(n);
% C = X(next,:)';

[C(:,1), index(1)] = datasample(X,1,1);

% 2) Compute distance from all other points to this centroid


% 3) choose one new data point at random as a new center using a weighted
% prob distribution where a point x is chosen with probability proportional
% to D(x)

minDist = ones(n,1);

% Select the rest of the seeds by a probabilistic model


for cen = 2:k
    % remove the chosen point
    X(index(cen-1),:) = []; maxSim(index(cen-1)) = [];
    
    minDist = min(mindist,1-(X*C(:,cen-1)).^2);
    
    sumDist = sum(maxSim);
    
    sampleProbability = minDist/sumDist;
    if ~isempty(stream)
        [C(:,cen), index(cen)] = datasample(stream,X,1,1,'Replace',false,...
            'Weights',sampleProbability);
    else
        [C(:,cen), index(cen)] = datasample(X,1,1,'Replace',false,...
            'Weights',sampleProbability);
    end
    
end



end
