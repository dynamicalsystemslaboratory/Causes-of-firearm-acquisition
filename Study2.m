clear

%load data
load BChecks% Background checks at a state level
load Mshooting_Patterns_states% Mass shooting at a state level (binary)
load laws%law restrictiveness index for each State
load population% population of each state


nboots=20000;% number of bootstrap replications

% state nº 7 is Connecticut and State nº 11 is Hawaii
population([7 11])=[]; % we eliminate Connecticut and Hawaii from population
BChecks(:,[7 11])=[];% we eliminate Connecticut and Hawaii from Bacground checks
laws([7 11])=[];% we eliminate Connecticut and Hawaii from law index

Mshooting_Patterns_states(:,[7 11])=[];% we eliminate Connecticut and Hawaii from mass shooting

% We use k-means to calculate 2 clusters according to law index, namely
% permissive and restrictive.
[idx,C] = kmeans(laws,2,'Distance','sqeuclidean');
a1=find(idx==1);
a2=find(idx==2);
if length(a1)<length(a2)
indres=a1;
indper=a2;
else
indres=a2;
indper=a1;
end;

popures=population(indres);%population of the restrictive states
popuper=population(indper);%population of the permissive states

sBChecks=BChecks>0;% we transform to binary Bacground checks
sBCres=sBChecks(:,indres);% Bacground checks for restrictive states (binary)
sBCper=sBChecks(:,indper);% Bacground checks for permissive states (binary)

sMS=Mshooting_Patterns_states;%rename the matrix Mshooting_Patterns_states into sMS
sMSres=Mshooting_Patterns_states(:,indres);% Mass shooting for restrictive states
sMSper=Mshooting_Patterns_states(:,indper);% Mass shooting for permissive states
sBChecks=BChecks>0;% we transform to binary Bacground checks



%%%%% Analysis 3 nodes

%Computation of transfer entropies for restrictive states TE_{X->Y_i|Z}
for s=1:length(indres)
    ms=sMS;
    ms(:,indres(s))=[];
    sms=sum(ms,2)>0;
    sms(140)=1;%include mass shooting conneticut 2010-08
    sms(168)=1;%include mass shooting conneticut 2012-12
    sms(11)=1;%include mass shooting Hawaii 1999-11
    sms(177)=1;% include mass shooting WDC 2013-09

    TEMStoBCrcondsms(s,1)=guns_ent2([sBCres(2:end,s) sms(1:end-1) sBCres(1:end-1,s)])-...
         guns_ent2([sms(1:end-1) sBCres(1:end-1,s)])-...
         guns_ent2([sBCres(2:end,s) sMSres(1:end-1,s) sms(1:end-1) sBCres(1:end-1,s)])+...
         guns_ent2([sMSres(1:end-1,s) sms(1:end-1) sBCres(1:end-1,s)]);
     
    TEsmstoBCrcondMS(s,1)=guns_ent2([sBCres(2:end,s) sMSres(1:end-1,s) sBCres(1:end-1,s)])-...
         guns_ent2([sMSres(1:end-1,s) sBCres(1:end-1,s)])-...
         guns_ent2([sBCres(2:end,s) sMSres(1:end-1,s) sms(1:end-1) sBCres(1:end-1,s)])+...
         guns_ent2([sMSres(1:end-1,s) sms(1:end-1) sBCres(1:end-1,s)]);
end;

%Computation of transfer entropies for permissive states

for s=1:length(indper)
    ms=sMS;
    ms(:,indper(s))=[];
    sms=sum(ms,2)>0;
    sms(140)=1;%include mass shooting conneticut 2010-08
    sms(168)=1;%include mass shooting conneticut 2012-12
    sms(11)=1;%include mass shooting Hawaii 1999-11
    sms(177)=1;% include mass shooting WDC 2013-09
    
    TEMStoBCpcondsms(s,1)=guns_ent2([sBCper(2:end,s) sms(1:end-1) sBCper(1:end-1,s)])-...
         guns_ent2([sms(1:end-1) sBCper(1:end-1,s)])-...
         guns_ent2([sBCper(2:end,s) sMSper(1:end-1,s) sms(1:end-1) sBCper(1:end-1,s)])+...
         guns_ent2([sMSper(1:end-1,s) sms(1:end-1) sBCper(1:end-1,s)]);
    TEsmstoBCpcondMS(s,1)=guns_ent2([sBCper(2:end,s) sMSper(1:end-1,s) sBCper(1:end-1,s)])-...
         guns_ent2([sMSper(1:end-1,s) sBCper(1:end-1,s)])-...
         guns_ent2([sBCper(2:end,s) sMSper(1:end-1,s) sms(1:end-1) sBCper(1:end-1,s)])+...
         guns_ent2([sMSper(1:end-1,s) sms(1:end-1) sBCper(1:end-1,s)]);
end;

%Computation of TE_{X->Y|Z} for the restrictive states
me_TEMStoBCrcondsms=1/sum(popures.^2)*sum(popures.*sqrt(TEMStoBCrcondsms))^2;
me_TEsmstoBCrcondMS=1/sum(popures.^2)*sum(popures.*sqrt(TEsmstoBCrcondMS))^2;

%Computation of TE_{X->Y|Z} for the permissive states
me_TEMStoBCpcondsms=1/sum(popuper.^2)*sum(popuper.*sqrt(TEMStoBCpcondsms))^2;
me_TEsmstoBCpcondMS=1/sum(popuper.^2)*sum(popuper.*sqrt(TEsmstoBCpcondMS))^2;

%%%%% Statistical analysis. Surrogate distribution

%Computation of the surrogate distribution of TE_{X->Y_i|Z} for permissive
%states
for s=1:length(indper)
    ms=sMS;
    ms(:,indper(s))=[];
    sms=sum(ms,2)>0;
    sms(140)=1;%include mass shooting conneticut 2010-08
    sms(168)=1;%include mass shooting conneticut 2012-12
    sms(11)=1;%include mass shooting Hawaii 1999-11
    sms(177)=1;% include mass shooting WDC 2013-09
    
    vector=[sMSper(:,s) sms sBCper(:,s)];
        parfor b=1:nboots
             
        MSb=guns_boot([vector(:,1) vector(:,2) vector(:,3)]);
TEMStoBCpcondsms(b,1)=guns_ent2([vector(2:end,3) vector(1:end-1,2) vector(1:end-1,3)])-...
         guns_ent2([vector(1:end-1,2) vector(1:end-1,3)])-...
         guns_ent2([vector(2:end,3) MSb(1:end-1) vector(1:end-1,2) vector(1:end-1,3)])+...
         guns_ent2([MSb(1:end-1) vector(1:end-1,2) vector(1:end-1,3)]);
     
        smsb=guns_boot([vector(:,2) vector(:,1) vector(:,3)]);
TEsmstoBCpcondMS(b,1)=guns_ent2([vector(2:end,3) vector(1:end-1,1) vector(1:end-1,3)])-...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,3)])-...
         guns_ent2([vector(2:end,3) vector(1:end-1,1) smsb(1:end-1) vector(1:end-1,3)])+...
         guns_ent2([vector(1:end-1,1) smsb(1:end-1) vector(1:end-1,3)]);
        end;
        
        TEMStoBCpcondsmsb{s}=TEMStoBCpcondsms;
        TEsmstoBCpcondMSb{s}=TEsmstoBCpcondMS;
        
 end;

 %Computation of the surrogate distribution of TE_{X->Y|Z} for permissive
 %states
 mTEMStoBCpcondsmsb=0;
 mTEsmstoBCpcondMSb=0;
 for  s=1:length(indper)
mTEMStoBCpcondsmsb=mTEMStoBCpcondsmsb+popuper(s)/sqrt(sum(popuper.^2))*sqrt(TEMStoBCpcondsmsb{s});
mTEsmstoBCpcondMSb=mTEsmstoBCpcondMSb+popuper(s)/sqrt(sum(popuper.^2))*sqrt(TEsmstoBCpcondMSb{s});
 end;

 %Calculate the percentile 95 of the surrogate distribution for permissive states
 cutTEMStoBCpcondsmsb=prctile(mTEMStoBCpcondsmsb.^2,95);
 cutTEsmstoBCpcondMSb=prctile(mTEsmstoBCpcondMSb.^2,95);

 %Calculate the p-value for significance of transfer entropy for permissive states
 pvalTEMStoBCpcondsmsb=1-sum(mTEMStoBCpcondsmsb.^2<me_TEMStoBCpcondsms)/nboots;
 pvalTEsmstoBCpcondMSb=1-sum(mTEsmstoBCpcondMSb.^2<me_TEsmstoBCpcondMS)/nboots;

 %Computation of the surrogate distribution of TE_{X->Y_i|Z} for restrictive states
 
for s=1:length(indres)
    ms=sMS;
    ms(:,indres(s))=[];
    sms=sum(ms,2)>0;
    sms(140)=1;%include mass shooting conneticut 2010-08
    sms(168)=1;%include mass shooting conneticut 2012-12
    sms(11)=1;%include mass shooting Hawaii 1999-11
    sms(177)=1;% include mass shooting WDC 2013-09
    
    vector=[sMSres(:,s) sms sBCres(:,s)];
      parfor b=1:nboots
            
         MSb=guns_boot([vector(:,1) vector(:,2) vector(:,3)]);
TEMStoBCrcondsms(b,1)=guns_ent2([vector(2:end,3) vector(1:end-1,2) vector(1:end-1,3)])-...
         guns_ent2([vector(1:end-1,2) vector(1:end-1,3)])-...
         guns_ent2([vector(2:end,3) MSb(1:end-1) vector(1:end-1,2) vector(1:end-1,3)])+...
         guns_ent2([MSb(1:end-1) vector(1:end-1,2) vector(1:end-1,3)]);
     
        smsb=guns_boot([vector(:,2) vector(:,1) vector(:,3)]);
TEsmstoBCrcondMS(b,1)=guns_ent2([vector(2:end,3) vector(1:end-1,1) vector(1:end-1,3)])-...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,3)])-...
         guns_ent2([vector(2:end,3) vector(1:end-1,1) smsb(1:end-1) vector(1:end-1,3)])+...
         guns_ent2([vector(1:end-1,1) smsb(1:end-1) vector(1:end-1,3)]);
        end;
        TEMStoBCrcondsmsb{s}=TEMStoBCrcondsms;
        TEsmstoBCrcondMSb{s}=TEsmstoBCrcondMS;
 end;
 
 %Computation of the surrogate distribution of TE_{X->Y|Z} for restrictive
%states
 mTEMStoBCrcondsmsb=0;
 mTEsmstoBCrcondMSb=0;
for  s=1:length(indres)
mTEMStoBCrcondsmsb=mTEMStoBCrcondsmsb+popures(s)/sqrt(sum(popures.^2))*sqrt(TEMStoBCrcondsmsb{s});
mTEsmstoBCrcondMSb=mTEsmstoBCrcondMSb+popures(s)/sqrt(sum(popures.^2))*sqrt(TEsmstoBCrcondMSb{s});
end;

 %Calculate the percentile 95 of the surrogate distribution for restrictive
 %states
 cutTEMStoBCrcondsmsb=prctile(mTEMStoBCrcondsmsb.^2,95);
 cutTEsmstoBCrcondMSb=prctile(mTEsmstoBCrcondMSb.^2,95);

 
 %Calculate the p-value for significance of transfer entropy for restrictive
 %states
 pvalTEMStoBCrcondsmsb=1-sum(mTEMStoBCrcondsmsb.^2<me_TEMStoBCrcondsms)/nboots;
 pvalTEsmstoBCrcondMSb=1-sum(mTEsmstoBCrcondMSb.^2<me_TEsmstoBCrcondMS)/nboots;

 % we gather all results in resultsStu2
 resultsStu2=[me_TEMStoBCrcondsms me_TEMStoBCpcondsms me_TEsmstoBCrcondMS  me_TEsmstoBCpcondMS;...
     cutTEMStoBCrcondsmsb cutTEMStoBCpcondsmsb cutTEsmstoBCrcondMSb cutTEsmstoBCpcondMSb;...
     pvalTEMStoBCrcondsmsb pvalTEMStoBCpcondsmsb pvalTEsmstoBCrcondMSb  pvalTEsmstoBCpcondMSb];
 
 save resultsStu2 resultsStu2
     
 