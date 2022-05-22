function H = entropy(X)
% EnTROPY   Estimate the entropy using
%           a cool binless estimator.
%
% Note [d,n] = size(X) (d dims, n samples).
%
% J Victor (2002) "Binless strategies for
% estimation of information from neural
% data", Physical Review E

% Use the machine precision for neighbor
% calculations.
precision = eps;
    
[d, n] = size(X);
% Calculate the nearest neighbor for each point.
for i=1:n
    % Remove the i'th vector.
    X_temp = X(:, find([1:n] - i));
        
    % Subtract off the i'th vector from all others.
    X_diff = X_temp - repmat(X(:, i), 1, n-1);
    
    % Calculate the minimum Euclidean distance.
    lambda(i) = min(sqrt(sum((X_diff).^2, 1)));
    
    % Ensure the distance is not zero.
    if (lambda(i) < precision)
        lambda(i) = precision;
    end
end
    
% The "Euler-Mascheroni" constant.
em = 0.5772156649015;
    
% Calculate area the of an d-dimensional sphere.
area = d * pi^(d/2) / gamma(1 + d/2);
    
% Calculate an estimate of entropy based on the 
% mean nearest neighbor distance using an equation
% from the above citation.
K = log2(area) + log2((n-1) / d) + em / log(2);
H = d * mean(log2(lambda)) + K;
    
   