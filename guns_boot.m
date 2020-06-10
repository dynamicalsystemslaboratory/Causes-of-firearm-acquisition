function results=guns_boot(x)
%the first column of x is the causal variable (X)
% need function permn
[T,K]=size(x);

symb=permn([0 1],K-1);%define all possible symbols
xc=x(:,2:K);%the matrix formed by the effect variable (Y) and conditioning variables (Z)
xx=x(:,1);% the effect variable (X)

for i=1:2^(K-1)
    loc=ismember(xc,symb(i,:),'rows');%binary vector with a 1 in the position where 
    %the symbol "symb(i,:)" appears in the row of xc
    pos{i}=find(loc>0);% keep the positions of the 2^{K-1} symbols
    clear loc   
end;

xxs=zeros(T,1);
for i=1:2^(K-1)
    xs=xx(pos{i});% this is the portion to be shuffled within the locations of the same symbol
    xxs(pos{i})=xs(randperm(length(pos{i})));%we fill in the shuffle causal vector
end;

results=xxs;%Bootstrap causal vector for the computation of the surrogate transfer entropy distribution.






