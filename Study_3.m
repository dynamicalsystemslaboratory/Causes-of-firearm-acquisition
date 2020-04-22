clear

% Load data
load BChecks % Background checks
load MOfc % Media output on firearm laws and regulations, namely MOfc
load MOs % Media output on shootings excluding firearm laws and regulations, namely MOs
load Mshooting_Patterns % Mass shooting
load D % Matrix of distances between the centroid of the states
load laws % Law restrictiveness index for each State
load population % Population of each State

nboots=20000; % Number of bootstrap replications

% State nº 7 is Connecticut and State nº 11 is Hawaii
population([7 11])=[]; % We eliminate Connecticut and Hawaii from population
BChecks(:,[7 11])=[]; % We eliminate Connecticut and Hawaii from Background checks
laws([7 11])=[]; % We eliminate Connecticut and Hawaii from law index
D(:,[7 11])=[]; % We eliminate Connecticut and Hawaii from the columns of the distance matrix
D([7 11],:)=[]; % We eliminate Connecticut and Hawaii from the rows of the distance matrix

MS=Mshooting_Patterns>0; % We transform to binary mass shooting
sBChecks=BChecks>0; % We transform to binary Bacground checks
sMOs=MOs>median(MOs); % We transform to binary MOs
sMOfc=MOfc>median(MOfc); % We transform to binary MOfc

% We use k-means to calculate 2 clusters according to law index, namely
% permissive and restrictive.
idx = kmeans(laws,2,'Distance','sqeuclidean');
a1=find(idx==1);
a2=find(idx==2);
% We know from inspection of the clusters that restrictive States
% are less than permissive.
% We need to include the following as initial center positions are random
if length(a1)<length(a2)
    indres=a1;
    indper=a2;
else
    indres=a2;
    indper=a1;
end

popures=population(indres); % Population of the restrictive States
popuper=population(indper); % Population of the permissive States

sBCres=sBChecks(:,indres); % Background checks for restrictive States (binary)
sBCper=sBChecks(:,indper); % Background checks for permissive States (binary)


nsizes=[1 3 5 7 9]; % Vector of different numbers of neighbors to be considered

for v=1:length(nsizes)
    n=nsizes(v)
    
    % Compute the connectivity matrix of n-nearest neighbors based on the
    % distance matrix D
    Wgeo=zeros(48,48);
    for i=1:48
        sdi=sort(D(i,:));
        indi = NaN(n+1,1); % n+1 as we pick up the State itself as the closest neighbor
        for j=1:n+1
            indi(j)=find(D(i,:)==sdi(j));
        end
        Wgeo(i,indi)=1;
    end
    Wgeo=Wgeo-eye(48); % Take out the State itself as a neighbor
    
    
    % Calculate the most frequent value of background checks (binary) in the n nearest neighbors
    % for the restrictive States
    geosBCres=zeros(228,length(indres));
    for t=1:228
        for i=1:length(indres)
            geosBCres(t,i)=mode(sBChecks(t,find(Wgeo(indres(i),:)>0)));
        end
    end
    % Calculate the most frequent value of background checks (binary) in the n nearest neighbors
    % for the permissive States
    geosBCper=zeros(228,length(indper));
    for t=1:228
        for i=1:length(indper)
            geosBCper(t,i)=mode(sBChecks(t,find(Wgeo(indper(i),:)>0)));
        end
    end
    
    
    
    
    
    %%%%% Analysis 3 nodes
    
    % Computation of transfer entropies for restrictive States TE_{X->Y_i|Z}
    TEBCvtoBCrcondMOfc = NaN(length(indres),1); % Preallocation
    TEMOfctoBCrcondBCv = NaN(length(indres),1); % Preallocation
    for s=1:length(indres)
        
        TEBCvtoBCrcondMOfc(s)=guns_ent2([sBCres(2:end,s) sBCres(1:end-1,s) sMOfc(1:end-1)])-...
            guns_ent2([sBCres(1:end-1,s) sMOfc(1:end-1)])-...
            guns_ent2([sBCres(2:end,s) sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOfc(1:end-1)])+...
            guns_ent2([sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOfc(1:end-1)]);
        
        TEMOfctoBCrcondBCv(s)=guns_ent2([sBCres(2:end,s) sBCres(1:end-1,s) geosBCres(1:end-1,s)])-...
            guns_ent2([sBCres(1:end-1,s) geosBCres(1:end-1,s)])-...
            guns_ent2([sBCres(2:end,s) sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOfc(1:end-1)])+...
            guns_ent2([sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOfc(1:end-1)]);
    end
    
    % Computation of transfer entropies for permissive States
    TEBCvtoBCpcondMOfc = NaN(length(indper),1); % Preallocation
    TEMOfctoBCpcondBCv = NaN(length(indper),1); % Preallocation
    for s=1:length(indper)
        
        TEBCvtoBCpcondMOfc(s)=guns_ent2([sBCper(2:end,s) sBCper(1:end-1,s) sMOfc(1:end-1)])-...
            guns_ent2([sBCper(1:end-1,s) sMOfc(1:end-1)])-...
            guns_ent2([sBCper(2:end,s) sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOfc(1:end-1)])+...
            guns_ent2([sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOfc(1:end-1)]);
        
        TEMOfctoBCpcondBCv(s)=guns_ent2([sBCper(2:end,s) sBCper(1:end-1,s) geosBCper(1:end-1,s)])-...
            guns_ent2([sBCper(1:end-1,s) geosBCper(1:end-1,s)])-...
            guns_ent2([sBCper(2:end,s) sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOfc(1:end-1)])+...
            guns_ent2([sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOfc(1:end-1)]);
    end
    
    % Computation of TE_{X->Y|Z} for the restrictive States
    me_TEBCvtoBCrcondMOfc=1/sum(popures.^2)*sum(popures.*sqrt(TEBCvtoBCrcondMOfc))^2;
    me_TEMOfctoBCrcondBCv=1/sum(popures.^2)*sum(popures.*sqrt(TEMOfctoBCrcondBCv))^2;
    % Computation of TE_{X->Y|Z} for the permissive States
    me_TEBCvtoBCpcondMOfc=1/sum(popuper.^2)*sum(popuper.*sqrt(TEBCvtoBCpcondMOfc))^2;
    me_TEMOfctoBCpcondBCv=1/sum(popuper.^2)*sum(popuper.*sqrt(TEMOfctoBCpcondBCv))^2;
    
    %%%%% Statistical analysis. Surrogate distribution
    
    % Computation of the surrogate distribution of TE_{X->Y_i|Z} for permissive
    % States
    
    for s=1:length(indper)
        
        vector=[sBCper(:,s) geosBCper(:,s) sMOfc];
        TEBCvtoBCpcondMOfc = NaN(nboots,1); % Preallocation
        TEMOfctoBCpcondBCv = NaN(nboots,1); % Preallocation
        parfor b=1:nboots
            
            geosBCperb=guns_boot([vector(:,2) vector(:,1) vector(:,3)]);
            TEBCvtoBCpcondMOfc(b)=guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,3)])-...
                guns_ent2([vector(1:end-1,1) vector(1:end-1,3)])-...
                guns_ent2([vector(2:end,1) vector(1:end-1,1) geosBCperb(1:end-1) vector(1:end-1,3)])+...
                guns_ent2([vector(1:end-1,1) geosBCperb(1:end-1) vector(1:end-1,3)]);
            
            sMOfcb=guns_boot([vector(:,3) vector(:,1) vector(:,2)]);
            TEMOfctoBCpcondBCv(b)=guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,2)])-...
                guns_ent2([vector(1:end-1,1) vector(1:end-1,2)])-...
                guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,2) sMOfcb(1:end-1)])+...
                guns_ent2([vector(1:end-1,1) vector(1:end-1,2) sMOfcb(1:end-1)]);
        end
        
        TEBCvtoBCpcondMOfcb{s}=TEBCvtoBCpcondMOfc;
        TEMOfctoBCpcondBCvb{s}=TEMOfctoBCpcondBCv;
        
    end
    
    % Computation of the surrogate distribution of TE_{X->Y|Z} for permissive
    % States
    mTEBCvtoBCpcondMOfcb=0;
    mTEMOfctoBCpcondBCvb=0;
    for  s=1:length(indper)
        mTEBCvtoBCpcondMOfcb=mTEBCvtoBCpcondMOfcb+popuper(s)/sqrt(sum(popuper.^2))*sqrt(TEBCvtoBCpcondMOfcb{s});
        mTEMOfctoBCpcondBCvb=mTEMOfctoBCpcondBCvb+popuper(s)/sqrt(sum(popuper.^2))*sqrt(TEMOfctoBCpcondBCvb{s});
    end
    
    % Calculate the percentile 95 of the surrogate distribution for permissive
    % States
    cutTEBCvtoBCpcondMOfcb=prctile(mTEBCvtoBCpcondMOfcb.^2,95);
    cutTEMOfctoBCpcondBCvb=prctile(mTEMOfctoBCpcondBCvb.^2,95);
    
    % Calculate the p-value for significance of transfer entropy for permissive
    % States
    pvalTEBCvtoBCpcondMOfcb=1-sum(mTEBCvtoBCpcondMOfcb.^2<me_TEBCvtoBCpcondMOfc)/nboots;
    pvalTEMOfctoBCpcondBCvb=1-sum(mTEMOfctoBCpcondBCvb.^2<me_TEMOfctoBCpcondBCv)/nboots;
    
    % Computation of the surrogate distribution of TE_{X->Y_i|Z} for restrictive states
    for s=1:length(indres)
        
        vector=[sBCres(:,s) geosBCres(:,s) sMOfc];
        TEBCvtoBCrcondMOfc = NaN(nboots,1); % Preallocation
        TEMOfctoBCrcondBCv = NaN(nboots,1); % Preallocation
        parfor b=1:nboots
            
            geosBCresb=guns_boot([vector(:,2) vector(:,1) vector(:,3)]);
            TEBCvtoBCrcondMOfc(b)=guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,3)])-...
                guns_ent2([vector(1:end-1,1) vector(1:end-1,3)])-...
                guns_ent2([vector(2:end,1) vector(1:end-1,1) geosBCresb(1:end-1) vector(1:end-1,3)])+...
                guns_ent2([vector(1:end-1,1) geosBCresb(1:end-1) vector(1:end-1,3)]);
            
            sMOfcb=guns_boot([vector(:,3) vector(:,1) vector(:,2)]);
            TEMOfctoBCrcondBCv(b)=guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,2)])-...
                guns_ent2([vector(1:end-1,1) vector(1:end-1,2)])-...
                guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,2) sMOfcb(1:end-1)])+...
                guns_ent2([vector(1:end-1,1) vector(1:end-1,2) sMOfcb(1:end-1)]);
        end
        TEBCvtoBCrcondMOfcb{s}=TEBCvtoBCrcondMOfc;
        TEMOfctoBCrcondBCvb{s}=TEMOfctoBCrcondBCv;
        
    end
    
    % Computation of the surrogate distribution of TE_{X->Y|Z} for restrictive
    % States
    mTEBCvtoBCrcondMOfcb=0;
    mTEMOfctoBCrcondBCvb=0;
    for s=1:length(indres)
        mTEBCvtoBCrcondMOfcb=mTEBCvtoBCrcondMOfcb+popures(s)/sqrt(sum(popures.^2))*sqrt(TEBCvtoBCrcondMOfcb{s});
        mTEMOfctoBCrcondBCvb=mTEMOfctoBCrcondBCvb+popures(s)/sqrt(sum(popures.^2))*sqrt(TEMOfctoBCrcondBCvb{s});
    end
    
    % Calculate the percentile 95 of the surrogate distribution for restrictive
    % States
    cutTEBCvtoBCrcondMOfcb=prctile(mTEBCvtoBCrcondMOfcb.^2,95);
    cutTEMOfctoBCrcondBCvb=prctile(mTEMOfctoBCrcondBCvb.^2,95);
    
    % Calculate the p-value for significance of transfer entropy for restrictive
    % States
    pvalTEBCvtoBCrcondMOfcb=1-sum(mTEBCvtoBCrcondMOfcb.^2<me_TEBCvtoBCrcondMOfc)/nboots;
    pvalTEMOfctoBCrcondBCvb=1-sum(mTEMOfctoBCrcondBCvb.^2<me_TEMOfctoBCrcondBCv)/nboots;
    
    
    results_all_n{v}=[me_TEBCvtoBCrcondMOfc me_TEBCvtoBCpcondMOfc me_TEMOfctoBCrcondBCv me_TEMOfctoBCpcondBCv;...
        cutTEBCvtoBCrcondMOfcb cutTEBCvtoBCpcondMOfcb  cutTEMOfctoBCrcondBCvb cutTEMOfctoBCpcondBCvb;...
        pvalTEBCvtoBCrcondMOfcb pvalTEBCvtoBCpcondMOfcb pvalTEMOfctoBCrcondBCvb pvalTEMOfctoBCpcondBCvb];
    
    clear me_TEBCvtoBCrcondMOfc
    clear me_TEBCvtoBCpcondMOfc
    clear me_TEMOfctoBCrcondBCv
    clear me_TEMOfctoBCpcondBCv
    clear TEBCvtoBCrcondMOfc
    clear TEBCvtoBCpcondMOfc
    clear TEMOfctoBCrcondBCv
    clear TEMOfctoBCpcondBCv
end
