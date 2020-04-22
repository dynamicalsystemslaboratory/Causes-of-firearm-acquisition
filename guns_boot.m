function results=guns_boot(x)
% The first column of x is the causal variable (X)
% Need function permn: https://www.mathworks.com/matlabcentral/fileexchange/7147-permn
[T,K]=size(x);

symb=permn([0 1],K-1); % Define all possible symbols for effect variable (Y) and conditioning variables (Z)
xc=x(:,2:K); % The matrix formed by the effect variable (Y) and conditioning variables (Z)
xx=x(:,1); % The causal variable (X)

for i=1:2^(K-1)
    loc=ismember(xc,symb(i,:),'rows'); % Binary vector with a 1 in the position where 
    % the symbol "symb(i,:)" appears in the row of xc
    pos{i}=find(loc>0); % Keep the positions of the 2^{K-1} symbols
    clear loc   
end

xxs=NaN(T,1);
for i=1:2^(K-1)
    xs=xx(pos{i}); % This is the portion to be shuffled within the locations of the same symbol
    xxs(pos{i})=xs(randperm(length(pos{i}))); % We fill in the shuffle causal vector
end

results=xxs; % Bootstrap causal vector for the computation of the surrogate transfer entropy distribution.






