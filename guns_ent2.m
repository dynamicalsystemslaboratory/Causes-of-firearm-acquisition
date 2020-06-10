function results=guns_ent2(x)
%Computation of entropy of the binary matrix X
%x is a T x K matrix

T=length(x);
[l,r,s]=unique(x,'rows'); %determine the vector s of different 2^K symbols in the rows os x
 pxs=hist(s,1:length(l))/T;% Compute the probability distribution of the 2^K symbols
 pxs=pxs+(pxs==0);% transforms 0 to 1
 logpxs=log2(pxs);% compute the log of the probabilities
 results=-sum(pxs.*logpxs);%compute the Shannon entropy

