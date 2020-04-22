clear

% Load data
load BChecks % Background checks
load MOfc % Media output on firearm laws and regulations, namely MOfc
load MOs % Media output on shootings excluding firearm laws and regulations, namely MOs
load Mshooting_Patterns % Mass shooting
load laws % Law restrictiveness index for each State
load population % Population of each State

nboots=20000; % Number of bootstrap replications

% State nº 7 is Connecticut and State nº 11 is Hawaii
population([7 11])=[]; % We eliminate Connecticut and Hawaii from population
BChecks(:,[7 11])=[]; % We eliminate Connecticut and Hawaii from Background checks
laws([7 11])=[]; % We eliminate Connecticut and Hawaii from law index

MS=Mshooting_Patterns>0; % We transform to binary mass shooting
sBChecks=BChecks>0; % We transform to binary Background checks
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




%%%%% Analysis 4 nodes

% Computation of transfer entropies for restrictive States TE_{X->Y_i|Z}
TEMStoBCrcondMOfc_MOs = NaN(length(indres),1); % Preallocation
TEMOfctoBCrcondMS_MOs = NaN(length(indres),1); % Preallocation
TEMOstoBCrcondMS_MOfc = NaN(length(indres),1); % Preallocation
for s=1:length(indres)
    
    TEMStoBCrcondMOfc_MOs(s)=guns_ent2([sBCres(2:end,s) sBCres(1:end-1,s) sMOfc(1:end-1) sMOs(1:end-1)])-...
        guns_ent2([sBCres(1:end-1,s) sMOfc(1:end-1) sMOs(1:end-1)])-...
        guns_ent2([sBCres(2:end,s) MS(1:end-1) sBCres(1:end-1,s) sMOfc(1:end-1) sMOs(1:end-1)])+...
        guns_ent2([MS(1:end-1) sBCres(1:end-1,s) sMOfc(1:end-1) sMOs(1:end-1)]);
    
    TEMOfctoBCrcondMS_MOs(s)=guns_ent2([sBCres(2:end,s) MS(1:end-1) sBCres(1:end-1,s) sMOs(1:end-1)])-...
        guns_ent2([MS(1:end-1) sBCres(1:end-1,s) sMOs(1:end-1)])-...
        guns_ent2([sBCres(2:end,s) MS(1:end-1) sBCres(1:end-1,s) sMOfc(1:end-1) sMOs(1:end-1)])+...
        guns_ent2([MS(1:end-1) sBCres(1:end-1,s) sMOfc(1:end-1) sMOs(1:end-1)]);
    
    TEMOstoBCrcondMS_MOfc(s)=guns_ent2([sBCres(2:end,s)  MS(1:end-1) sBCres(1:end-1,s) sMOfc(1:end-1)])-...
        guns_ent2([MS(1:end-1) sBCres(1:end-1,s) sMOfc(1:end-1)])-...
        guns_ent2([sBCres(2:end,s) MS(1:end-1) sBCres(1:end-1,s) sMOfc(1:end-1) sMOs(1:end-1)])+...
        guns_ent2([MS(1:end-1) sBCres(1:end-1,s) sMOfc(1:end-1) sMOs(1:end-1)]);
end

% Computation of transfer entropies for permissive States
TEMStoBCpcondMOfc_MOs = NaN(length(indper),1); % Preallocation
TEMOfctoBCpcondMS_MOs = NaN(length(indper),1); % Preallocation
TEMOstoBCpcondMS_MOfc = NaN(length(indper),1); % Preallocation
for s=1:length(indper)
    
    TEMStoBCpcondMOfc_MOs(s)=guns_ent2([sBCper(2:end,s) sBCper(1:end-1,s) sMOfc(1:end-1) sMOs(1:end-1)])-...
        guns_ent2([sBCper(1:end-1,s) sMOfc(1:end-1) sMOs(1:end-1)])-...
        guns_ent2([sBCper(2:end,s) MS(1:end-1) sBCper(1:end-1,s) sMOfc(1:end-1) sMOs(1:end-1)])+...
        guns_ent2([MS(1:end-1) sBCper(1:end-1,s) sMOfc(1:end-1) sMOs(1:end-1)]);
    
    TEMOfctoBCpcondMS_MOs(s)=guns_ent2([sBCper(2:end,s)  MS(1:end-1) sBCper(1:end-1,s) sMOs(1:end-1)])-...
        guns_ent2([MS(1:end-1) sBCper(1:end-1,s) sMOs(1:end-1)])-...
        guns_ent2([sBCper(2:end,s) MS(1:end-1) sBCper(1:end-1,s) sMOfc(1:end-1) sMOs(1:end-1)])+...
        guns_ent2([MS(1:end-1) sBCper(1:end-1,s) sMOfc(1:end-1) sMOs(1:end-1)]);
    
    TEMOstoBCpcondMS_MOfc(s)=guns_ent2([sBCper(2:end,s) MS(1:end-1) sBCper(1:end-1,s) sMOfc(1:end-1)])-...
        guns_ent2([MS(1:end-1) sBCper(1:end-1,s) sMOfc(1:end-1)])-...
        guns_ent2([sBCper(2:end,s) MS(1:end-1) sBCper(1:end-1,s) sMOfc(1:end-1) sMOs(1:end-1)])+...
        guns_ent2([MS(1:end-1) sBCper(1:end-1,s) sMOfc(1:end-1) sMOs(1:end-1)]);
end

% Computation of TE_{X->Y|Z} for the restrictive States
me_TEMStoBCrcondMOfc_MOs=1/sum(popures.^2)*sum(popures.*sqrt(TEMStoBCrcondMOfc_MOs))^2;
me_TEMOfctoBCrcondMS_MOs=1/sum(popures.^2)*sum(popures.*sqrt(TEMOfctoBCrcondMS_MOs))^2;
me_TEMOstoBCrcondMS_MOfc=1/sum(popures.^2)*sum(popures.*sqrt(TEMOstoBCrcondMS_MOfc))^2;

% Computation of TE_{X->Y|Z} for the permissive States
me_TEMStoBCpcondMOfc_MOs=1/sum(popuper.^2)*sum(popuper.*sqrt(TEMStoBCpcondMOfc_MOs))^2;
me_TEMOfctoBCpcondMS_MOs=1/sum(popuper.^2)*sum(popuper.*sqrt(TEMOfctoBCpcondMS_MOs))^2;
me_TEMOstoBCpcondMS_MOfc=1/sum(popuper.^2)*sum(popuper.*sqrt(TEMOstoBCpcondMS_MOfc))^2;

%%%%% Statistical analysis. Surrogate distribution

% Computation of the surrogate distribution of TE_{X->Y_i|Z} for permissive
% States
for s=1:length(indper)
    vector=[MS sBCper(:,s) sMOfc sMOs];
    TEMStoBCpcondMOfc_MOs = NaN(nboots,1); % Preallocation
    TEMOfctoBCpcondMS_MOs = NaN(nboots,1); % Preallocation
    TEMOstoBCpcondMS_MOfc = NaN(nboots,1); % Preallocation
    parfor b=1:nboots
        
        MSb=guns_boot([vector(:,1) vector(:,2) vector(:,3) vector(:,4)]);
        TEMStoBCpcondMOfc_MOs(b)=guns_ent2([vector(2:end,2) vector(1:end-1,2) vector(1:end-1,3) vector(1:end-1,4)])-...
            guns_ent2([vector(1:end-1,2) vector(1:end-1,3) vector(1:end-1,4)])-...
            guns_ent2([vector(2:end,2) MSb(1:end-1) vector(1:end-1,2) vector(1:end-1,3) vector(1:end-1,4)])+...
            guns_ent2([MSb(1:end-1) vector(1:end-1,2) vector(1:end-1,3) vector(1:end-1,4)]);
        
        sMOfcb=guns_boot([vector(:,3) vector(:,1) vector(:,2) vector(:,4)]);
        TEMOfctoBCpcondMS_MOs(b)=guns_ent2([vector(2:end,2)  vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,4)])-...
            guns_ent2([vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,4)])-...
            guns_ent2([vector(2:end,2) vector(1:end-1,1) vector(1:end-1,2) sMOfcb(1:end-1) vector(1:end-1,4)])+...
            guns_ent2([vector(1:end-1,1) vector(1:end-1,2) sMOfcb(1:end-1) vector(1:end-1,4)]);
        
        sMOsb=guns_boot([vector(:,4) vector(:,1) vector(:,2) vector(:,3)]);
        TEMOstoBCpcondMS_MOfc(b)=guns_ent2([vector(2:end,2) vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,3)])-...
            guns_ent2([vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,3)])-...
            guns_ent2([vector(2:end,2) vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,3) sMOsb(1:end-1)])+...
            guns_ent2([vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,3) sMOsb(1:end-1)]);
    end
    TEMStoBCpcondMOfc_MOsb{s}=TEMStoBCpcondMOfc_MOs;
    TEMOfctoBCpcondMS_MOsb{s}=TEMOfctoBCpcondMS_MOs;
    TEMOstoBCpcondMS_MOfcb{s}=TEMOstoBCpcondMS_MOfc;
end

% Computation of the surrogate distribution of TE_{X->Y|Z} for permissive
% States
mTEMStoBCpcondMOfc_MOsb=0;
mTEMOfctoBCpcondMS_MOsb=0;
mTEMOstoBCpcondMS_MOfcb=0;
for s=1:length(indper)
    mTEMStoBCpcondMOfc_MOsb=mTEMStoBCpcondMOfc_MOsb+popuper(s)/sqrt(sum(popuper.^2))*sqrt(TEMStoBCpcondMOfc_MOsb{s});
    mTEMOfctoBCpcondMS_MOsb=mTEMOfctoBCpcondMS_MOsb+popuper(s)/sqrt(sum(popuper.^2))*sqrt(TEMOfctoBCpcondMS_MOsb{s});
    mTEMOstoBCpcondMS_MOfcb=mTEMOstoBCpcondMS_MOfcb+popuper(s)/sqrt(sum(popuper.^2))*sqrt(TEMOstoBCpcondMS_MOfcb{s});
end

% Calculate the percentile 95 of the surrogate distribution for permissive
% States
cutTEMStoBCpcondMOfc_MOsb=prctile(mTEMStoBCpcondMOfc_MOsb.^2,95);
cutTEMOfctoBCpcondMS_MOsb=prctile(mTEMOfctoBCpcondMS_MOsb.^2,95);
cutTEMOstoBCpcondMS_MOfcb=prctile(mTEMOstoBCpcondMS_MOfcb.^2,95);

% Calculate the p-value for significance of transfer entropy for permissive
% States
pvalTEMStoBCpcondMOfc_MOsb=1-sum(mTEMStoBCpcondMOfc_MOsb.^2<me_TEMStoBCpcondMOfc_MOs)/nboots;
pvalTEMOfctoBCpcondMS_MOsb=1-sum(mTEMOfctoBCpcondMS_MOsb.^2<me_TEMOfctoBCpcondMS_MOs)/nboots;
pvalTEMOstoBCpcondMS_MOfcb=1-sum(mTEMOstoBCpcondMS_MOfcb.^2<me_TEMOstoBCpcondMS_MOfc)/nboots;

% Computation of the surrogate distribution of TE_{X->Y_i|Z} for restrictive States

for s=1:length(indres)
    vector=[MS sBCres(:,s) sMOfc sMOs];
    TEMStoBCrcondMOfc_MOs = NaN(nboots,1); % Preallocation
    TEMOfctoBCrcondMS_MOs = NaN(nboots,1); % Preallocation
    TEMOstoBCrcondMS_MOfc = NaN(nboots,1); % Preallocation
    parfor b=1:nboots
        
        MSb=guns_boot([vector(:,1) vector(:,2) vector(:,3) vector(:,4)]);
        TEMStoBCrcondMOfc_MOs(b)=guns_ent2([vector(2:end,2) vector(1:end-1,2) vector(1:end-1,3) vector(1:end-1,4)])-...
            guns_ent2([vector(1:end-1,2) vector(1:end-1,3) vector(1:end-1,4)])-...
            guns_ent2([vector(2:end,2) MSb(1:end-1) vector(1:end-1,2) vector(1:end-1,3) vector(1:end-1,4)])+...
            guns_ent2([MSb(1:end-1) vector(1:end-1,2) vector(1:end-1,3) vector(1:end-1,4)]);
        
        sMOfcb=guns_boot([vector(:,3) vector(:,1) vector(:,2) vector(:,4)]);
        TEMOfctoBCrcondMS_MOs(b)=guns_ent2([vector(2:end,2)  vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,4)])-...
            guns_ent2([vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,4)])-...
            guns_ent2([vector(2:end,2) vector(1:end-1,1) vector(1:end-1,2) sMOfcb(1:end-1) vector(1:end-1,4)])+...
            guns_ent2([vector(1:end-1,1) vector(1:end-1,2) sMOfcb(1:end-1) vector(1:end-1,4)]);
        
        sMOsb=guns_boot([vector(:,4) vector(:,1) vector(:,2) vector(:,3)]);
        TEMOstoBCrcondMS_MOfc(b)=guns_ent2([vector(2:end,2) vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,3)])-...
            guns_ent2([vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,3)])-...
            guns_ent2([vector(2:end,2) vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,3) sMOsb(1:end-1)])+...
            guns_ent2([vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,3) sMOsb(1:end-1)]);
    end
    TEMStoBCrcondMOfc_MOsb{s}=TEMStoBCrcondMOfc_MOs;
    TEMOfctoBCrcondMS_MOsb{s}=TEMOfctoBCrcondMS_MOs;
    TEMOstoBCrcondMS_MOfcb{s}=TEMOstoBCrcondMS_MOfc;
end

% Computation of the surrogate distribution of TE_{X->Y|Z} for restrictive
% States
mTEMStoBCrcondMOfc_MOsb=0;
mTEMOfctoBCrcondMS_MOsb=0;
mTEMOstoBCrcondMS_MOfcb=0;
for s=1:length(indres)
    mTEMStoBCrcondMOfc_MOsb=mTEMStoBCrcondMOfc_MOsb+popures(s)/sqrt(sum(popures.^2))*sqrt(TEMStoBCrcondMOfc_MOsb{s});
    mTEMOfctoBCrcondMS_MOsb=mTEMOfctoBCrcondMS_MOsb+popures(s)/sqrt(sum(popures.^2))*sqrt(TEMOfctoBCrcondMS_MOsb{s});
    mTEMOstoBCrcondMS_MOfcb=mTEMOstoBCrcondMS_MOfcb+popures(s)/sqrt(sum(popures.^2))*sqrt(TEMOstoBCrcondMS_MOfcb{s});
end

% Calculate the percentile 95 of the surrogate distribution for restrictive
% States
cutTEMStoBCrcondMOfc_MOsb=prctile(mTEMStoBCrcondMOfc_MOsb.^2,95);
cutTEMOfctoBCrcondMS_MOsb=prctile(mTEMOfctoBCrcondMS_MOsb.^2,95);
cutTEMOstoBCrcondMS_MOfcb=prctile(mTEMOstoBCrcondMS_MOfcb.^2,95);

% Calculate the p-value for significance of transfer entropy for restrictive
% States
pvalTEMStoBCrcondMOfc_MOsb=1-sum(mTEMStoBCrcondMOfc_MOsb.^2<me_TEMStoBCrcondMOfc_MOs)/nboots;
pvalTEMOfctoBCrcondMS_MOsb=1-sum(mTEMOfctoBCrcondMS_MOsb.^2<me_TEMOfctoBCrcondMS_MOs)/nboots;
pvalTEMOstoBCrcondMS_MOfcb=1-sum(mTEMOstoBCrcondMS_MOfcb.^2<me_TEMOstoBCrcondMS_MOfc)/nboots;