function results=guns_ent2(x)
% Computation of entropy of the binary matrix X
% x is a T x K matrix

T=length(x);
[l,r,s]=unique(x,'rows'); % Determine the vector s of different 2^K symbols in the rows os x
pxs=hist(s,1:length(l))/T; % Compute the probability distribution of the 2^K symbols
pxs=pxs+(pxs==0); % Transforms 0 to 1 to avoid -infinity when taking the logarithm
logpxs=log2(pxs); % Compute the log of the probabilities
results=-sum(pxs.*logpxs); % Compute the Shannon entropy

