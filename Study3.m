clear

%load data
load BChecks% Background checks
load MOs_Patterns%media output on shootings excluding firearm laws and regulations, namely MOs
load MOfc_Patterns%media output on firearm laws and regulations, namely MOfc
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

sBChecks=BChecks>0;% we transform to binary Bacground checks

sMOs=mode(serMOs,2);%binary time series for MOs
sMOfc=mode(serMOfc,2);%binary time series of MOfc

nboots=1000;% number of bootstrap replications

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

    TEBCvtoBCrcondMOs_MOfc(s,1)=guns_ent2([sBCres(2:end,s) sBCres(1:end-1,s) sMOs(1:end-1) sMOfc(1:end-1)])-...
         guns_ent2([sBCres(1:end-1,s) sMOs(1:end-1) sMOfc(1:end-1)])-...
         guns_ent2([sBCres(2:end,s) sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOs(1:end-1) sMOfc(1:end-1)])+...
         guns_ent2([sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOs(1:end-1) sMOfc(1:end-1)]);
     
    TEMOstoBCrcondBCv_MOfc(s,1)=guns_ent2([sBCres(2:end,s) sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOfc(1:end-1)])-...
         guns_ent2([sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOfc(1:end-1)])-...
         guns_ent2([sBCres(2:end,s) sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOs(1:end-1) sMOfc(1:end-1)])+...
         guns_ent2([sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOs(1:end-1) sMOfc(1:end-1)]);
     
    TEMOfctoBCrcondBCv_MOs(s,1)=guns_ent2([sBCres(2:end,s) sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOs(1:end-1)])-...
         guns_ent2([sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOs(1:end-1)])-...
         guns_ent2([sBCres(2:end,s) sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOs(1:end-1) sMOfc(1:end-1)])+...
         guns_ent2([sBCres(1:end-1,s) geosBCres(1:end-1,s) sMOs(1:end-1) sMOfc(1:end-1)]);
end;

%Computation of transfer entropies for permissive states
for s=1:length(indper)

    TEBCvtoBCpcondMOs_MOfc(s,1)=guns_ent2([sBCper(2:end,s) sBCper(1:end-1,s) sMOs(1:end-1) sMOfc(1:end-1)])-...
         guns_ent2([sBCper(1:end-1,s) sMOs(1:end-1) sMOfc(1:end-1)])-...
         guns_ent2([sBCper(2:end,s) sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOs(1:end-1) sMOfc(1:end-1)])+...
         guns_ent2([sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOs(1:end-1) sMOfc(1:end-1)]);
     
    TEMOstoBCpcondBCv_MOfc(s,1)=guns_ent2([sBCper(2:end,s) sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOfc(1:end-1)])-...
         guns_ent2([sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOfc(1:end-1)])-...
         guns_ent2([sBCper(2:end,s) sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOs(1:end-1) sMOfc(1:end-1)])+...
         guns_ent2([sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOs(1:end-1) sMOfc(1:end-1)]);
     
    TEMOfctoBCpcondBCv_MOs(s,1)=guns_ent2([sBCper(2:end,s) sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOs(1:end-1)])-...
         guns_ent2([sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOs(1:end-1)])-...
         guns_ent2([sBCper(2:end,s) sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOs(1:end-1) sMOfc(1:end-1)])+...
         guns_ent2([sBCper(1:end-1,s) geosBCper(1:end-1,s) sMOs(1:end-1) sMOfc(1:end-1)]);
end;


%Computation of TE_{X->Y|Z} for the restrictive states
  
me_TEBCvtoBCrcondMOs_MOfc=1/sum(popures.^2)*sum(popures.*sqrt(TEBCvtoBCrcondMOs_MOfc))^2;
me_TEMOstoBCrcondBCv_MOfc=1/sum(popures.^2)*sum(popures.*sqrt(TEMOstoBCrcondBCv_MOfc))^2;
me_TEMOfctoBCrcondBCv_MOs=1/sum(popures.^2)*sum(popures.*sqrt(TEMOfctoBCrcondBCv_MOs))^2;
%Computation of TE_{X->Y|Z} for the permissive states
  
me_TEBCvtoBCpcondMOs_MOfc=1/sum(popuper.^2)*sum(popuper.*sqrt(TEBCvtoBCpcondMOs_MOfc))^2;
me_TEMOstoBCpcondBCv_MOfc=1/sum(popuper.^2)*sum(popuper.*sqrt(TEMOstoBCpcondBCv_MOfc))^2;
me_TEMOfctoBCpcondBCv_MOs=1/sum(popuper.^2)*sum(popuper.*sqrt(TEMOfctoBCpcondBCv_MOs))^2;
%%%%% Statistical analysis. Surrogate distribution

%Computation of the surrogate distribution of TE_{X->Y_i|Z} for permissive
%states

for s=1:length(indper)

    vector=[sBCper(:,s) geosBCper(:,s) sMOs sMOfc];
        parfor b=1:nboots
             
        geosBCperb=guns_boot([vector(:,2) vector(:,1) vector(:,3) vector(:,4)]);
TEBCvtoBCpcondMOs_MOfc(b,1)=guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,3) vector(1:end-1,4)])-...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,3) vector(1:end-1,4)])-...
         guns_ent2([vector(2:end,1) vector(1:end-1,1) geosBCperb(1:end-1) vector(1:end-1,3) vector(1:end-1,4)])+...
         guns_ent2([vector(1:end-1,1) geosBCperb(1:end-1) vector(1:end-1,3) vector(1:end-1,4)]);
     
        sMOsb=guns_boot([vector(:,3) vector(:,1) vector(:,2) vector(:,4)]);
TEMOstoBCpcondBCv_MOfc(b,1)=guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,4)])-...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,4)])-...
         guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,2) sMOsb(1:end-1) vector(1:end-1,4)])+...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,2) sMOsb(1:end-1) vector(1:end-1,4)]);
             sMOfcb=guns_boot([vector(:,3) vector(:,1) vector(:,2) vector(:,4)]);
          
        sMOfcb=guns_boot([vector(:,4) vector(:,1) vector(:,2) vector(:,3)]);
TEMOfctoBCpcondBCv_MOs(b,1)=guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,3)])-...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,3)])-...
         guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,3) sMOfcb(1:end-1)])+...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,2)  vector(1:end-1,3) sMOfcb(1:end-1)]);
        end;
          
        TEBCvtoBCpcondMOs_MOfcb{s}=TEBCvtoBCpcondMOs_MOfc;
        TEMOstoBCpcondBCv_MOfcb{s}=TEMOstoBCpcondBCv_MOfc;
        TEMOfctoBCpcondBCv_MOsb{s}=TEMOfctoBCpcondBCv_MOs;
        
 end;

 %Computation of the surrogate distribution of TE_{X->Y|Z} for permissive
 %states
 mTEBCvtoBCpcondMOs_MOfcb=0;
 mTEMOstoBCpcondBCv_MOfcb=0;
 mTEMOfctoBCpcondBCv_MOsb=0;
 for  s=1:length(indper)
mTEBCvtoBCpcondMOs_MOfcb=mTEBCvtoBCpcondMOs_MOfcb+popuper(s)/sqrt(sum(popuper.^2))*sqrt(TEBCvtoBCpcondMOs_MOfcb{s});
mTEMOstoBCpcondBCv_MOfcb=mTEMOstoBCpcondBCv_MOfcb+popuper(s)/sqrt(sum(popuper.^2))*sqrt(TEMOstoBCpcondBCv_MOfcb{s});
mTEMOfctoBCpcondBCv_MOsb=mTEMOfctoBCpcondBCv_MOsb+popuper(s)/sqrt(sum(popuper.^2))*sqrt(TEMOfctoBCpcondBCv_MOsb{s});

 end;
 
 
 %Calculate the percentile 95 of the surrogate distribution for permissive
 %states
 cutTEBCvtoBCpcondMOs_MOfcb=prctile(mTEBCvtoBCpcondMOs_MOfcb.^2,95);
 cutTEMOstoBCpcondBCv_MOfcb=prctile(mTEMOstoBCpcondBCv_MOfcb.^2,95);
 cutTEMOfctoBCpcondBCv_MOsb=prctile(mTEMOfctoBCpcondBCv_MOsb.^2,95);

 %  
 %Calculate the p-value for significance of transfer entropy for permissive
 %states
 pvalTEBCvtoBCpcondMOs_MOfcb=1-sum(mTEBCvtoBCpcondMOs_MOfcb.^2<me_TEBCvtoBCpcondMOs_MOfc)/nboots;
 pvalTEMOstoBCpcondBCv_MOfcb=1-sum(mTEMOstoBCpcondBCv_MOfcb.^2<me_TEMOstoBCpcondBCv_MOfc)/nboots;
 pvalTEMOfctoBCpcondBCv_MOsb=1-sum(mTEMOfctoBCpcondBCv_MOsb.^2<me_TEMOfctoBCpcondBCv_MOs)/nboots;

 %Computation of the surrogate distribution of TE_{X->Y_i|Z} for restrictive states
for s=1:length(indres)

    vector=[sBCres(:,s) geosBCres(:,s) sMOs sMOfc];
        parfor b=1:nboots
             
        geosBCresb=guns_boot([vector(:,2) vector(:,1) vector(:,3) vector(:,4)]);
TEBCvtoBCrcondMOs_MOfc(b,1)=guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,3) vector(1:end-1,4)])-...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,3) vector(1:end-1,4)])-...
         guns_ent2([vector(2:end,1) vector(1:end-1,1) geosBCresb(1:end-1) vector(1:end-1,3) vector(1:end-1,4)])+...
         guns_ent2([vector(1:end-1,1) geosBCresb(1:end-1) vector(1:end-1,3) vector(1:end-1,4)]);
     
        sMOsb=guns_boot([vector(:,3) vector(:,1) vector(:,2) vector(:,4)]);
TEMOstoBCrcondBCv_MOfc(b,1)=guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,4)])-...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,4)])-...
         guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,2) sMOsb(1:end-1) vector(1:end-1,4)])+...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,2) sMOsb(1:end-1) vector(1:end-1,4)]);
             sMOfcb=guns_boot([vector(:,3) vector(:,1) vector(:,2) vector(:,4)]);
          
        sMOfcb=guns_boot([vector(:,4) vector(:,1) vector(:,2) vector(:,3)]);
TEMOfctoBCrcondBCv_MOs(b,1)=guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,3)])-...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,3)])-...
         guns_ent2([vector(2:end,1) vector(1:end-1,1) vector(1:end-1,2) vector(1:end-1,3) sMOfcb(1:end-1)])+...
         guns_ent2([vector(1:end-1,1) vector(1:end-1,2)  vector(1:end-1,3) sMOfcb(1:end-1)]);
        end;
          
        TEBCvtoBCrcondMOs_MOfcb{s}=TEBCvtoBCrcondMOs_MOfc;
        TEMOstoBCrcondBCv_MOfcb{s}=TEMOstoBCrcondBCv_MOfc;
        TEMOfctoBCrcondBCv_MOsb{s}=TEMOfctoBCrcondBCv_MOs;
        
 end;

 %Computation of the surrogate distribution of TE_{X->Y|Z} for resmissive
 %states
 mTEBCvtoBCrcondMOs_MOfcb=0;
 mTEMOstoBCrcondBCv_MOfcb=0;
 mTEMOfctoBCrcondBCv_MOsb=0;
 for  s=1:length(indres)
mTEBCvtoBCrcondMOs_MOfcb=mTEBCvtoBCrcondMOs_MOfcb+popures(s)/sqrt(sum(popures.^2))*sqrt(TEBCvtoBCrcondMOs_MOfcb{s});
mTEMOstoBCrcondBCv_MOfcb=mTEMOstoBCrcondBCv_MOfcb+popures(s)/sqrt(sum(popures.^2))*sqrt(TEMOstoBCrcondBCv_MOfcb{s});
mTEMOfctoBCrcondBCv_MOsb=mTEMOfctoBCrcondBCv_MOsb+popures(s)/sqrt(sum(popures.^2))*sqrt(TEMOfctoBCrcondBCv_MOsb{s});

 end;
 
 
 %Calculate the rescentile 95 of the surrogate distribution for resmissive
 %states
 cutTEBCvtoBCrcondMOs_MOfcb=prctile(mTEBCvtoBCrcondMOs_MOfcb.^2,95);
 cutTEMOstoBCrcondBCv_MOfcb=prctile(mTEMOstoBCrcondBCv_MOfcb.^2,95);
 cutTEMOfctoBCrcondBCv_MOsb=prctile(mTEMOfctoBCrcondBCv_MOsb.^2,95);

 
 %Calculate the p-value for significance of transfer entropy for resmissive
 %states
 pvalTEBCvtoBCrcondMOs_MOfcb=1-sum(mTEBCvtoBCrcondMOs_MOfcb.^2<me_TEBCvtoBCrcondMOs_MOfc)/nboots;
 pvalTEMOstoBCrcondBCv_MOfcb=1-sum(mTEMOstoBCrcondBCv_MOfcb.^2<me_TEMOstoBCrcondBCv_MOfc)/nboots;
 pvalTEMOfctoBCrcondBCv_MOsb=1-sum(mTEMOfctoBCrcondBCv_MOsb.^2<me_TEMOfctoBCrcondBCv_MOs)/nboots;

 
 results_all_n{v}=[me_TEBCvtoBCrcondMOs_MOfc me_TEBCvtoBCpcondMOs_MOfc me_TEMOstoBCrcondBCv_MOfc...
     me_TEMOstoBCpcondBCv_MOfc me_TEMOfctoBCrcondBCv_MOs me_TEMOfctoBCpcondBCv_MOs;...
     cutTEBCvtoBCrcondMOs_MOfcb cutTEBCvtoBCpcondMOs_MOfcb  cutTEMOstoBCrcondBCv_MOfcb...
     cutTEMOstoBCpcondBCv_MOfcb cutTEMOfctoBCrcondBCv_MOsb cutTEMOfctoBCpcondBCv_MOsb;...
     pvalTEBCvtoBCrcondMOs_MOfcb pvalTEBCvtoBCpcondMOs_MOfcb pvalTEMOstoBCrcondBCv_MOfcb...
     pvalTEMOstoBCpcondBCv_MOfcb pvalTEMOfctoBCrcondBCv_MOsb pvalTEMOfctoBCpcondBCv_MOsb];
 
 save results_all_n results_all_n
  
     clear me_TEBCvtoBCrcondMOs_MOfc 
     clear me_TEBCvtoBCpcondMOs_MOfc 
     clear me_TEMOstoBCrcondBCv_MOfc
     clear me_TEMOstoBCpcondBCv_MOfc 
     clear me_TEMOfctoBCrcondBCv_MOs 
     clear me_TEMOfctoBCpcondBCv_MOs
     clear TEBCvtoBCrcondMOs_MOfc
	 clear TEBCvtoBCpcondMOs_MOfc
     clear TEMOstoBCrcondBCv_MOfc
     clear TEMOstoBCpcondBCv_MOfc 
     clear TEMOfctoBCrcondBCv_MOs 
     clear TEMOfctoBCpcondBCv_MOs
end;

 
