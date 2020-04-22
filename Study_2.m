clear

% Load data
load BChecks % Background checks at a State level
load Mshooting_Patterns_states % Mass shooting at a State level (binary)
load laws % Law restrictiveness index for each State
load population % population of each State


nboots=20000; % Number of bootstrap replications

% State nº 7 is Connecticut and State nº 11 is Hawaii
population([7 11])=[]; % We eliminate Connecticut and Hawaii from population
BChecks(:,[7 11])=[]; % We eliminate Connecticut and Hawaii from Baclground checks
laws([7 11])=[]; % We eliminate Connecticut and Hawaii from law index
Mshooting_Patterns_states(:,[7 11])=[]; % We eliminate Connecticut and Hawaii from mass shooting

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

sBChecks=BChecks>0; % We transform to binary Background checks
sBCres=sBChecks(:,indres); % Background checks for restrictive States (binary)
sBCper=sBChecks(:,indper); % Background checks for permissive States (binary)

sMS=Mshooting_Patterns_states; % Rename the matrix Mshooting_Patterns_states into sMS
sMSres=Mshooting_Patterns_states(:,indres); % Mass shooting for restrictive States
sMSper=Mshooting_Patterns_states(:,indper); % Mass shooting for permissive States


%%%%% Analysis 3 nodes

% Computation of transfer entropies for restrictive States TE_{X->Y_i|Z}
TEMStoBCrcondsms = NaN(length(indres),1); % Preallocation
TEsmstoBCrcondMS = NaN(length(indres),1); % Preallocation
for s=1:length(indres)
    ms=sMS;
    ms(:,indres(s))=[];
    sms=sum(ms,2)>0; % Mass shootings in other States
    
    TEMStoBCrcondsms(s)=guns_ent2([sBCres(2:end,s) sms(1:end-1) sBCres(1:end-1,s)])-...
        guns_ent2([sms(1:end-1) sBCres(1:end-1,s)])-...
        guns_ent2([sBCres(2:end,s) sMSres(1:end-1,s) sms(1:end-1) sBCres(1:end-1,s)])+...
        guns_ent2([sMSres(1:end-1,s) sms(1:end-1) sBCres(1:end-1,s)]);
    
    TEsmstoBCrcondMS(s)=guns_ent2([sBCres(2:end,s) sMSres(1:end-1,s) sBCres(1:end-1,s)])-...
        guns_ent2([sMSres(1:end-1,s) sBCres(1:end-1,s)])-...
        guns_ent2([sBCres(2:end,s) sMSres(1:end-1,s) sms(1:end-1) sBCres(1:end-1,s)])+...
        guns_ent2([sMSres(1:end-1,s) sms(1:end-1) sBCres(1:end-1,s)]);
end

% Computation of transfer entropies for permissive States
TEMStoBCpcondsms = NaN(length(indper),1); % Preallocation
TEsmstoBCpcondMS = NaN(length(indper),1); % Preallocation
for s=1:length(indper)
    ms=sMS;
    ms(:,indper(s))=[];
    sms=sum(ms,2)>0;
    
    TEMStoBCpcondsms(s)=guns_ent2([sBCper(2:end,s) sms(1:end-1) sBCper(1:end-1,s)])-...
        guns_ent2([sms(1:end-1) sBCper(1:end-1,s)])-...
        guns_ent2([sBCper(2:end,s) sMSper(1:end-1,s) sms(1:end-1) sBCper(1:end-1,s)])+...
        guns_ent2([sMSper(1:end-1,s) sms(1:end-1) sBCper(1:end-1,s)]);
    
    TEsmstoBCpcondMS(s)=guns_ent2([sBCper(2:end,s) sMSper(1:end-1,s) sBCper(1:end-1,s)])-...
        guns_ent2([sMSper(1:end-1,s) sBCper(1:end-1,s)])-...
        guns_ent2([sBCper(2:end,s) sMSper(1:end-1,s) sms(1:end-1) sBCper(1:end-1,s)])+...
        guns_ent2([sMSper(1:end-1,s) sms(1:end-1) sBCper(1:end-1,s)]);
end

% Computation of TE_{X->Y|Z} for the restrictive States
me_TEMStoBCrcondsms=1/sum(popures.^2)*sum(popures.*sqrt(TEMStoBCrcondsms))^2;
me_TEsmstoBCrcondMS=1/sum(popures.^2)*sum(popures.*sqrt(TEsmstoBCrcondMS))^2;

% Computation of TE_{X->Y|Z} for the permissive States
me_TEMStoBCpcondsms=1/sum(popuper.^2)*sum(popuper.*sqrt(TEMStoBCpcondsms))^2;
me_TEsmstoBCpcondMS=1/sum(popuper.^2)*sum(popuper.*sqrt(TEsmstoBCpcondMS))^2;

%%%%% Statistical analysis. Surrogate distribution

% Computation of the surrogate distribution of TE_{X->Y_i|Z} for permissive
% States
for s=1:length(indper)
    ms=sMS;
    ms(:,indper(s))=[];
    sms=sum(ms,2)>0;
    vector=[sMSper(:,s) sms sBCper(:,s)];
    TEMStoBCpcondsms = NaN(nboots,1); % Preallocation
    TEsmstoBCpcondMS = NaN(nboots,1); % Preallocation
    parfor b=1:nboots
        
        MSb=guns_boot([vector(:,1) vector(:,2) vector(:,3)]);
        TEMStoBCpcondsms(b)=guns_ent2([vector(2:end,3) vector(1:end-1,2) vector(1:end-1,3)])-...
            guns_ent2([vector(1:end-1,2) vector(1:end-1,3)])-...
            guns_ent2([vector(2:end,3) MSb(1:end-1) vector(1:end-1,2) vector(1:end-1,3)])+...
            guns_ent2([MSb(1:end-1) vector(1:end-1,2) vector(1:end-1,3)]);
        
        smsb=guns_boot([vector(:,2) vector(:,1) vector(:,3)]);
        TEsmstoBCpcondMS(b)=guns_ent2([vector(2:end,3) vector(1:end-1,1) vector(1:end-1,3)])-...
            guns_ent2([vector(1:end-1,1) vector(1:end-1,3)])-...
            guns_ent2([vector(2:end,3) vector(1:end-1,1) smsb(1:end-1) vector(1:end-1,3)])+...
            guns_ent2([vector(1:end-1,1) smsb(1:end-1) vector(1:end-1,3)]);
    end
    
    TEMStoBCpcondsmsb{s}=TEMStoBCpcondsms;
    TEsmstoBCpcondMSb{s}=TEsmstoBCpcondMS;
    
end

% Computation of the surrogate distribution of TE_{X->Y|Z} for permissive
% States
mTEMStoBCpcondsmsb=0;
mTEsmstoBCpcondMSb=0;
for s=1:length(indper)
    mTEMStoBCpcondsmsb=mTEMStoBCpcondsmsb+popuper(s)/sqrt(sum(popuper.^2))*sqrt(TEMStoBCpcondsmsb{s});
    mTEsmstoBCpcondMSb=mTEsmstoBCpcondMSb+popuper(s)/sqrt(sum(popuper.^2))*sqrt(TEsmstoBCpcondMSb{s});
end

% Calculate the percentile 95 of the surrogate distribution for permissive States
cutTEMStoBCpcondsmsb=prctile(mTEMStoBCpcondsmsb.^2,95);
cutTEsmstoBCpcondMSb=prctile(mTEsmstoBCpcondMSb.^2,95);

% Calculate the p-value for significance of transfer entropy for permissive States
pvalTEMStoBCpcondsmsb=1-sum(mTEMStoBCpcondsmsb.^2<me_TEMStoBCpcondsms)/nboots;
pvalTEsmstoBCpcondMSb=1-sum(mTEsmstoBCpcondMSb.^2<me_TEsmstoBCpcondMS)/nboots;

% Computation of the surrogate distribution of TE_{X->Y_i|Z} for restrictive States

for s=1:length(indres)
    ms=sMS;
    ms(:,indres(s))=[];
    sms=sum(ms,2)>0;
    vector=[sMSres(:,s) sms sBCres(:,s)];
    TEMStoBCrcondsms = NaN(nboots,1); % Preallocation
    TEsmstoBCrcondMS = NaN(nboots,1); % Preallocation
    parfor b=1:nboots
        
        MSb=guns_boot([vector(:,1) vector(:,2) vector(:,3)]);
        TEMStoBCrcondsms(b)=guns_ent2([vector(2:end,3) vector(1:end-1,2) vector(1:end-1,3)])-...
            guns_ent2([vector(1:end-1,2) vector(1:end-1,3)])-...
            guns_ent2([vector(2:end,3) MSb(1:end-1) vector(1:end-1,2) vector(1:end-1,3)])+...
            guns_ent2([MSb(1:end-1) vector(1:end-1,2) vector(1:end-1,3)]);
        
        smsb=guns_boot([vector(:,2) vector(:,1) vector(:,3)]);
        TEsmstoBCrcondMS(b)=guns_ent2([vector(2:end,3) vector(1:end-1,1) vector(1:end-1,3)])-...
            guns_ent2([vector(1:end-1,1) vector(1:end-1,3)])-...
            guns_ent2([vector(2:end,3) vector(1:end-1,1) smsb(1:end-1) vector(1:end-1,3)])+...
            guns_ent2([vector(1:end-1,1) smsb(1:end-1) vector(1:end-1,3)]);
    end
    TEMStoBCrcondsmsb{s}=TEMStoBCrcondsms;
    TEsmstoBCrcondMSb{s}=TEsmstoBCrcondMS;
end

% Computation of the surrogate distribution of TE_{X->Y|Z} for restrictive
% States
mTEMStoBCrcondsmsb=0;
mTEsmstoBCrcondMSb=0;
for s=1:length(indres)
    mTEMStoBCrcondsmsb=mTEMStoBCrcondsmsb+popures(s)/sqrt(sum(popures.^2))*sqrt(TEMStoBCrcondsmsb{s});
    mTEsmstoBCrcondMSb=mTEsmstoBCrcondMSb+popures(s)/sqrt(sum(popures.^2))*sqrt(TEsmstoBCrcondMSb{s});
end

% Calculate the percentile 95 of the surrogate distribution for restrictive
% States
cutTEMStoBCrcondsmsb=prctile(mTEMStoBCrcondsmsb.^2,95);
cutTEsmstoBCrcondMSb=prctile(mTEsmstoBCrcondMSb.^2,95);

% Calculate the p-value for significance of transfer entropy for restrictive
% States
pvalTEMStoBCrcondsmsb=1-sum(mTEMStoBCrcondsmsb.^2<me_TEMStoBCrcondsms)/nboots;
pvalTEsmstoBCrcondMSb=1-sum(mTEsmstoBCrcondMSb.^2<me_TEsmstoBCrcondMS)/nboots;

