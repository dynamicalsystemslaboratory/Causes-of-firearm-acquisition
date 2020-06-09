clear

%load data
load BChecks% Background checks
load MOs_Patterns%media output on shootings excluding firearm laws and regulations, namely MOs
load MOfc_Patterns%media output on firearm laws and regulations, namely MOfc
%load Mshooting_Patterns% Mass shooting
load MStotal_Patterns% Mass shooting
load laws%law restrictiveness index for each State
load population% population of each state
load D%distance matrix between centroides of states


mediansMOs=median(MOs_Patterns);%compute the median of MOs
mediansMOfc=median(MOfc_Patterns);%compute the mean of MOfc
for i=1:5
    serMOs(:,i)=MOs_Patterns(:,i)>mediansMOs(i);% symbolization of MOs time series
    serMOfc(:,i)=MOfc_Patterns(:,i)>mediansMOfc(i);% symbolization of MOfc time series
end

nboots=20000;% number of bootstrap replications

% state nº 7 is Connecticut and State nº 11 is Hawaii
population([7 11])=[]; % we eliminate Connecticut and Hawaii from population
BChecks(:,[7 11])=[];% we eliminate Connecticut and Hawaii from Bacground checks
laws([7 11])=[];% we eliminate Connecticut and Hawaii from law index
D(:,[7 11])=[];% we eliminate Connecticut and Hawaii 
D([7 11],:)=[];% we eliminate Connecticut and Hawaii 

MS=MStotal_Patterns;% rename mass shootings
sBChecks=BChecks>0;% we transform to binary Bacground checks

sMOs=mode(serMOs,2);%binary time series for MOs
sMOfc=mode(serMOfc,2);%binary time series of MOfc

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

sBCres=sBChecks(:,indres);% Bacground checks for restrictive states (binary)
sBCper=sBChecks(:,indper);% Bacground checks for permissive states (binary)


nsizes=[1 3 5 7 9];% vector of different number of neighbors to be considered

for v=1:length(nsizes)
   n=nsizes(v)
   
   % Compute the conectivity matrix of n-nearest neighbors based on the
   % distance matrix D
Wgeo=zeros(48,48);
for i=1:48
    sdi=sort(D(i,:));
    for j=1:n+1
    indi(j)=find(D(i,:)==sdi(j));
    end;
    Wgeo(i,indi)=1;
end;
    Wgeo=Wgeo-eye(48);
    

% calculate the most frequent value of background checks (binary) in the n nearest neighbors
% for the restrictive states
geosBCres=zeros(228,length(indres));
for t=1:228
    for i=1:length(indres)
        geosBCres(t,i)=mode(sBChecks(t,find(Wgeo(indres(i),:)>0)));
    end;
end;
% calculate the most frequent value of background checks (binary) in the n nearest neighbors
% for the permissive states
geosBCper=zeros(228,length(indper));
for t=1:228
    for i=1:length(indper)
 geosBCper(t,i)=mode(sBChecks(t,find(Wgeo(indper(i),:)>0)));
    end;
end;


%%%%% Analysis 3 nodes

%Computation of transfer entropies for restrictive states TE_{X->Y_i|Z}
for s=1:length(indres)

    TEBCvtoBCrcondMOfc(s,1)=guns_ent2([sBCres(2:end,s) sBCres(1:end-1,s) sMOfc(1:end-1)])-...
         guns_ent2([sBCres(1:end-1,s) sMOfc(1:end-1)])-...
         guns_ent2([sBCres(2:end,s) sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOfc(1:end-1)])+...
         guns_ent2([sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOfc(1:end-1)]);
     
    TEMOfctoBCrcondBCv(s,1)=guns_ent2([sBCres(2:end,s) sBCres(1:end-1,s) geosBCres(1:end-1,s)])-...
         guns_ent2([sBCres(1:end-1,s) geosBCres(1:end-1,s)])-...
         guns_ent2([sBCres(2:end,s) sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOfc(1:end-1)])+...
         guns_ent2([sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOfc(1:end-1)]);
end;

%Computation of transfer entropies for permissive states
for s=1:length(indper)

    TEBCvtoBCpcondMOfc(s,1)=guns_ent2([sBCper(2:end,s) sBCper(1:end-1,s) sMOfc(1:end-1)])-...
         guns_ent2([sBCper(1:end-1,s) sMOfc(1:end-1)])-...
         guns_ent2([sBCper(2:end,s) sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOfc(1:end-1)])+...
         guns_ent2([sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOfc(1:end-1)]);
     
    TEMOfctoBCpcondBCv(s,1)=guns_ent2([sBCper(2:end,s) sBCper(1:end-1,s) geosBCper(1:end-1,s)])-...
         guns_ent2([sBCper(1:end-1,s) geosBCper(1:end-1,s)])-...
         guns_ent2([sBCper(2:end,s) sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOfc(1:end-1)])+...
         guns_ent2([sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOfc(1:end-1)]);
end;

%Computation of TE_{X->Y|Z} for the restrictive states
me_TEBCvtoBCrcondMOfc=1/sum(popures.^2)*sum(popures.*sqrt(TEBCvtoBCrcondMOfc))^2;
me_TEMOfctoBCrcondBCv=1/sum(popures.^2)*sum(popures.*sqrt(TEMOfctoBCrcondBCv))^2;
%Computation of TE_{X->Y|Z} for the permissive states
me_TEBCvtoBCpcondMOfc=1/sum(popuper.^2)*sum(popuper.*sqrt(TEBCvtoBCpcondMOfc))^2;
me_TEMOfctoBCpcondBCv=1/sum(popuper.^2)*sum(popuper.*sqrt(TEMOfctoBCpcondBCv))^2;

%%%%% Statistical analysis. Surrogate distribution

%Computation of the surrogate distribution of TE_{X->Y_i|Z} for permissive
%states

for s=1:length(indper)

    vector=[sBCper(:,s) geosBCper(:,s) sMOfc];
        parfor b=1:nboots
             
        geosBCperb=guns_boot([vector(:,2) vector(:,1) vector(:,3)]);
TEBCvtoBCpcondMOfc(b,1)=guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,3)])-...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,3)])-...
         guns_ent2([vector(2:end,1) vector(1:end-1,1) geosBCperb(1:end-1) vector(1:end-1,3)])+...
         guns_ent2([vector(1:end-1,1) geosBCperb(1:end-1) vector(1:end-1,3)]);
     
        sMOfcb=guns_boot([vector(:,3) vector(:,1) vector(:,2)]);
TEMOfctoBCpcondBCv(b,1)=guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,2)])-...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,2)])-...
         guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,2) sMOfcb(1:end-1)])+...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,2) sMOfcb(1:end-1)]);
        end;
        
        TEBCvtoBCpcondMOfcb{s}=TEBCvtoBCpcondMOfc;
        TEMOfctoBCpcondBCvb{s}=TEMOfctoBCpcondBCv;
        
 end;

 %Computation of the surrogate distribution of TE_{X->Y|Z} for permissive
 %states
 mTEBCvtoBCpcondMOfcb=0;
 mTEMOfctoBCpcondBCvb=0;
 for  s=1:length(indper)
mTEBCvtoBCpcondMOfcb=mTEBCvtoBCpcondMOfcb+popuper(s)/sqrt(sum(popuper.^2))*sqrt(TEBCvtoBCpcondMOfcb{s});
mTEMOfctoBCpcondBCvb=mTEMOfctoBCpcondBCvb+popuper(s)/sqrt(sum(popuper.^2))*sqrt(TEMOfctoBCpcondBCvb{s});
 end;
 
 %Calculate the percentile 95 of the surrogate distribution for permissive
 %states
 cutTEBCvtoBCpcondMOfcb=prctile(mTEBCvtoBCpcondMOfcb.^2,95);
 cutTEMOfctoBCpcondBCvb=prctile(mTEMOfctoBCpcondBCvb.^2,95);
  
 %Calculate the p-value for significance of transfer entropy for permissive
 %states
 pvalTEBCvtoBCpcondMOfcb=1-sum(mTEBCvtoBCpcondMOfcb.^2<me_TEBCvtoBCpcondMOfc)/nboots;
 pvalTEMOfctoBCpcondBCvb=1-sum(mTEMOfctoBCpcondBCvb.^2<me_TEMOfctoBCpcondBCv)/nboots;

 %Computation of the surrogate distribution of TE_{X->Y_i|Z} for restrictive states
for s=1:length(indres)
    
  vector=[sBCres(:,s) geosBCres(:,s) sMOfc];
      parfor b=1:nboots
            
                   
        geosBCresb=guns_boot([vector(:,2) vector(:,1) vector(:,3)]);
TEBCvtoBCrcondMOfc(b,1)=guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,3)])-...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,3)])-...
         guns_ent2([vector(2:end,1) vector(1:end-1,1) geosBCresb(1:end-1) vector(1:end-1,3)])+...
         guns_ent2([vector(1:end-1,1) geosBCresb(1:end-1) vector(1:end-1,3)]);
     
        sMOfcb=guns_boot([vector(:,3) vector(:,1) vector(:,2)]);
TEMOfctoBCrcondBCv(b,1)=guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,2)])-...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,2)])-...
         guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,2) sMOfcb(1:end-1)])+...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,2) sMOfcb(1:end-1)]);
        end;
        TEBCvtoBCrcondMOfcb{s}=TEBCvtoBCrcondMOfc;
        TEMOfctoBCrcondBCvb{s}=TEMOfctoBCrcondBCv;
        
 end;

%Computation of the surrogate distribution of TE_{X->Y|Z} for restrictive
%states
 mTEBCvtoBCrcondMOfcb=0;
 mTEMOfctoBCrcondBCvb=0;
 for  s=1:length(indres)
mTEBCvtoBCrcondMOfcb=mTEBCvtoBCrcondMOfcb+popures(s)/sqrt(sum(popures.^2))*sqrt(TEBCvtoBCrcondMOfcb{s});
mTEMOfctoBCrcondBCvb=mTEMOfctoBCrcondBCvb+popures(s)/sqrt(sum(popures.^2))*sqrt(TEMOfctoBCrcondBCvb{s});
 end;
 
  %Calculate the percentile 95 of the surrogate distribution for restrictive
 %states
 cutTEBCvtoBCrcondMOfcb=prctile(mTEBCvtoBCrcondMOfcb.^2,95);
 cutTEMOfctoBCrcondBCvb=prctile(mTEMOfctoBCrcondBCvb.^2,95);
 
%Calculate the p-value for significance of transfer entropy for restrictive
 %states
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
end;

 