%% Graph Theory CS & DB, Oct 2017, @ King's College
% Loading cleaned DWI connectomes (ACT-based MRTRIX) 

%% Load data
clear all 
NsubsC=19; % Number of subjects group 1 
NsubsD=20; % Number of subjects group 2 
NumRois=90; 
Ntot=NsubsC+NsubsD

okC=[1:2,4:14,16:21] % exclude C03 and C15 (hemorrhage & shunt)
okD=[1:19,21] % exlcude D20, xtreme movement

for n=1:NsubsC
    disp(['Loading subject ',num2str(okC(n))]);
    tic;
    % Charlottes matlab only handles unzipped niis for load_nii
    dataClean(:,:,okC(n))=load(sprintf('/Volumes/DriveChar/Fossa_Posteriors/DWI/MergedFiles2/normalised/tckwtsAverageResponse/C%02d_sift_connectome_AR.csv',okC(n))); % Directory in your HD
    toc;
end
dataCleanAd=dataClean(:,:,okC)

for n=1:NsubsD
    disp(['Loading subject ',okD(n)]);
    tic;
    % Charlottes matlab only handles unzipped niis for load_nii
    dataClean2(:,:,okD(n))=load(sprintf('/Volumes/DriveChar/Fossa_Posteriors/DWI/MergedFiles2/normalised/tckwtsAverageResponse/D%02d_sift_connectome_AR.csv',okD(n))); % Directory in your HD
    toc;
end
dataClean2Ad=dataClean2(:,:,okD)

% Data complete
Data_all=cat(3, dataCleanAd,dataClean2Ad)
group=[zeros(1,NsubsC), ones(1,NsubsD)]

%% Load mu-values
for n=1:NsubsC
    disp(['Loading subject ',num2str(okC(n))]);
    tic;
    % Charlottes matlab only handles unzipped niis for load_nii
    MU(:,okC(n))=load(sprintf('/Volumes/DriveChar/Fossa_Posteriors/DWI/MergedFiles2/normalised/tcks/C%02d_tckweights_mu',okC(n))); % Directory in your HD
    toc;
end
MU1=MU(:,okC)

for n=1:NsubsD
    disp(['Loading subject ',num2str(okD(n))]);
    tic;
    % Charlottes matlab only handles unzipped niis for load_nii
    MUad(:,okD(n))=load(sprintf('/Volumes/DriveChar/Fossa_Posteriors/DWI/MergedFiles2/normalised/tcks/D%02d_tckweights_mu',okD(n))); % Directory in your HD
    toc;
end
MU2=MUad(:,okD)

Data_all_mu=cat(2, MU1,MU2)

%% Metrics for each cost  = Cost-weighted metrics (based on raw networks)
% ESTIMATE DENSITY:
maxDens=zeros(Ntot,1);
for n=1:size(Data_all,3)
    maxDens(n)=density_und((Data_all(1:90,1:90,n).*~eye(90))>0); 
end
netCostVec=0.80:-0.01:0.01; % Cost range between minimum network density value and 
Nth=length(netCostVec);

weightsToUse=[1]; % I only used trackweights now
Nconns=length(weightsToUse);

N=Ntot; Naal=NumRois;
globEffW=zeros(N,Nth);
locEffW=zeros(N,Nth);
entrop=zeros(N,Nth);
netStrength=zeros(N,Nth);
pathLength=zeros(N,Nth);
clust=zeros(N,Nth,Nconns);
globEffBin=zeros(N,Nth);
locEffBin=zeros(N,Nth);
pathLengthBin=zeros(N,Nth);
clustBin=zeros(N,Nth);
avgDegree=zeros(N,Nth);

netCost=zeros(N,Nth);
connToTest=Data_all(1:Naal,1:Naal,:)

%% In case of normalized values
% globEff_norm=zeros(N,length(netCostVec));
% locEff_norm=zeros(N,length(netCostVec));
% totalEnergy=zeros(N,length(netCostVec));
% clustBin=zeros(N,length(netCostVec));
% pathLengthBin=zeros(N,length(netCostVec));
for s=1:N    
    c=0;  
    for netCostAux=netCostVec
        c=c+1;  
        connCostAux=connToTest(:,:,s).*~eye(Naal);
        connAuxRescaled=connCostAux./max(max(max(Data_all(1:Naal,1:Naal,:))));
        % calculate thresholded connectome for subj s, with threshold netCostAux
        connCost=threshold_proportional(connAuxRescaled,netCostAux); 
%         connFA=connFAall(:,:,s).*~eye(Naal) % USE THIS CONNECOME FOR FA METRICS
        % create mask for these values
        connCostMask=double(connCost>0); 
        disp(['Processing subject ',num2str(s),', conn ',num2str(1),'/',num2str(Nconns),', network cost ', num2str(netCostAux),'...']);
          tic;  
        %for connT=1:Nconns                   
            connAux=connToTest(:,:,s).*~eye(Naal); 
            % multiply the mask and the subject's connectome
            connAux2=connAuxRescaled.*connCostMask;        
            % normalize the thresholded connectome
            connAuxNorm=2.*connAux2./sum(sum(connAux2)); 
                % for subject s, and cost c, give  GT metrics of normalized
                % networks
%                 pathLength_norm(s,c)=charpath(distance_wei(1./connAuxNorm));
%                 clust_norm(s,c)=mean(clustering_coef_wu(connAuxNorm)); 
%                 globEff_norm(s,c)=efficiency_wei(connAuxNorm);
%                 nodEff_norm(:,s,c)=efficiency_wei(connAuxNorm,1);
%                 locEff_norm(s,c)=mean(nodEff_norm(:,s,c));
%                 totalEnergy(s,c)=sum(sum(connAux2));                 
%                 globEffW(s,c)=efficiency_wei(connAux2);
%                 nodEffW(:,s,c)=efficiency_wei(connAux2,1);
%                 locEffW(s,c)=mean(nodEffW(:,s,c));
%                 globEffBin(s,c)=efficiency_bin(double(connAux2>0));
%                 locEffBin(s,c)=mean(efficiency_bin(double(connAux2>0),1));                 
                nodStrength(:,s,c)=strengths_und(connAux2);
                nodStrength_norm(:,s,c)=strengths_und(connAuxNorm);
%                avgStrength(s,c)=mean(nodStrength(:,s,c));                 
%                 nodDegree(:,s,c)=degrees_und(connAux2);
%                 avgDegree(s,c)=mean(nodDegree(:,s,c));      
%                   pathLengthBin(s,c)=charpath(distance_bin(double(connAux2>0)));
%                   clustBin(s,c)=mean(clustering_coef_bu(double(connAux2>0)));
        %end
        toc
    end
end
save('GT_metricsSIFTadapted_costANDnormalized_justStrengths')

%% Visualize metric against cost
%load('GT_metricsDWI_costANDnormalized')
% Choose which GT metric you want to display against cost
aux2=sigmaBin;
figTitle=' Smallworldness Binary ';
figure, hold on;
aux=aux2(:,:,connT);
aux(isnan(aux))=0;
aux(isinf(aux))=0;    
Vars1=aux;
filtNetCost=1:length(netCostVec);
filt=1:N;
% Visualize the GT metric Vars1 on Y-axis against cost (or network density)
% for both groups
plot(netCostVec(filtNetCost),mean(Vars1(filt(group==0),filtNetCost),1),'Color',[0 .50 .30], 'Linewidth',3);
plot(netCostVec(filtNetCost),mean(Vars1(filt(group==1),filtNetCost),1),'Color',[0.0 .0 .8], 'Linewidth',3);
plot(netCostVec(filtNetCost),mean(Vars1(filt(group==0),filtNetCost),1)+std(Vars1(filt(group==0),filtNetCost),1),'--','Color',[0 .50 .30],'Linewidth',3);
plot(netCostVec(filtNetCost),mean(Vars1(filt(group==0),filtNetCost),1)-std(Vars1(filt(group==0),filtNetCost),1),'--','Color',[0 .50 .30], 'Linewidth',3);
plot(netCostVec(filtNetCost),mean(Vars1(filt(group==1),filtNetCost),1),'Color',[0.0 .0 .8], 'Linewidth',3);
plot(netCostVec(filtNetCost),mean(Vars1(filt(group==1),filtNetCost),1)+std(Vars1(filt(group==1),filtNetCost),1),'--','Color',[0.2 .20 .8], 'Linewidth',3);
plot(netCostVec(filtNetCost),mean(Vars1(filt(group==1),filtNetCost),1)-std(Vars1(filt(group==1),filtNetCost),1),'--','Color',[0.2 .20 .8], 'Linewidth',3);
auxTit=' % Fibers';
title([figTitle,auxTit],'FontSize',16); leg=legend('Patients','Controls'); set(leg,'FontSize',20);
xlabel('Network density','FontSize',18); ylabel(figTitle)
set(gca,'FontSize',18);
hold off;
% additional stars for significance on figure costfunction
clear pval pvalAux n;
hold on;
for n=1:size(Vars1,2)
    pvalAux=anova1(Vars1(:,n), group,'off');
    pval(n)=pvalAux(1)
        if pval(n)<0.001
            plot(netCostVec(n),mean(Vars1(filt(group==0),n),1)./1.2,'*k','LineWidth',3,'MarkerSize',5);
        elseif pval(n)<0.01
            plot(netCostVec(n),mean(Vars1(filt(group==0),n),1)./1.2,'^k','LineWidth',3,'MarkerSize',5);
        elseif pval(n)<0.05
            plot(netCostVec(n),mean(Vars1(filt(group==0),n),1)./1.2,'+k','LineWidth',3,'MarkerSize',5);
        end
end

%% Cost-integrated selection
% Create figure of energy against cost for the two groups, to see how much
% of the total energy (or edge weights) is explained by which network cost
figure, plot(netCostVec(end:-1:1),(mean(totalEnergy(group==0,end:-1:1))/mean(totalEnergy(group==0,1))),'r')
hold on; 
plot(netCostVec(end:-1:1),(mean(totalEnergy(group==1,end:-1:1))/mean(totalEnergy(group==1,1))),'b');

% For all subjects:
figure, hold on;
for n=1:39
    plot(netCostVec(end:-1:1),((totalEnergy(n,end:-1:1))/mean(totalEnergy(n,1))),'k')
end
hold off;
% Check for lowest subject which value of cost explains 95% (in this case
% .26 = value 53 of netCostVec
for s=1:Ntot
LimCostInt_GEw(s)=mean(globEffW(s,55:80))
LimCostInt_LEw(s)=mean(locEffW(s,55:80))
LimCostInt_Netstrength(s)=mean(netStrength(s,55:80))
end
add=anova1(LimCostInt_GEw, group, 'Off');
bh = boxplot(LimCostInt_GEw, group, 'Labels',{'Patients','Controls'}); box 'off'
set(bh(:,:),'linewidth',2)
title('Cost-integrated Global Efficiency by group', 'FontSize',12); 
ylabel('Global Efficiency','FontSize',20); xlabel('Group','FontSize',18);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 18)
txt1 = ['p = ' num2str(add)] ;
txt2 = [' * '] 
text(1.4,max(LimCostInt_GEw),txt2, 'FontSize',30)

add=anova1(LimCostInt_LEw, group, 'Off');
bh = boxplot(LimCostInt_LEw, group, 'Labels',{'Patients','Controls'}); box 'off'
set(bh(:,:),'linewidth',2)
title('Cost-integrated Local Efficiency by group', 'FontSize',12); 
ylabel('Global Efficiency','FontSize',20); xlabel('Group','FontSize',18);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 18)
txt1 = ['p = ' num2str(add)] ;
txt2 = [' * '] 
text(1.4,max(LimCostInt_LEw),txt2, 'FontSize',30)

add=anova1(LimCostInt_Netstrength, group, 'Off');
bh = boxplot(LimCostInt_Netstrength, group, 'Labels',{'Patients','Controls'}); box 'off'
set(bh(:,:),'linewidth',2)
title('Cost-integrated Average Strength by group', 'FontSize',12); 
ylabel('Global Efficiency','FontSize',20); xlabel('Group','FontSize',18);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 18)
txt1 = ['p = ' num2str(add)] ;
txt2 = [' * '] 
text(1.4,max(LimCostInt_Netstrength),txt2, 'FontSize',30)

%% In case of using smallworldness, you need Clustering/clust random and path/path random
%load('GT_metricsSIFTadapted_costANDnormalized')
Nrand=10;
globEffR=zeros(N,Nth,Nrand,Nconns);
locEffR=zeros(N,Nth,Nrand,Nconns);
pathLengthR=zeros(N,Nth,Nrand,Nconns);
clustR=zeros(N,Nth,Nrand,Nconns);
globEffBinR=zeros(N,Nth,Nrand);
locEffBinR=zeros(N,Nth,Nrand);
pathLengthBinR=zeros(N,Nth,Nrand);
clustBinR=zeros(N,Nth,Nrand);

for s=1:N    
    c=0;
    for netCostAux=netCostVec
        tic;
        c=c+1;  
        connCostAux=connToTest(:,:,s).*~eye(Naal);
        % calculate thresholded connectome for subj s, with threshold netCostAux
        connCost=threshold_proportional(connCostAux,netCostAux); 
        % create mask for these values
        connCostMask=double(connCost>0); 
        for connT=1:Nconns                   
            disp(['Processing subject ',num2str(s),', conn ',num2str(connT),'/',num2str(Nconns),', network cost ', num2str(netCostAux),'...']);
            connAux=connToTest(:,:,s,connT).*~eye(Naal); 
            % multiply the mask and the subject's connectome
            connAux2=connAux.*connCostMask;        
            % normalize the thresholded connectome
            connAuxNorm=2.*connAux2./sum(sum(connAux2)); 
               for r=1:Nrand
                    % create randomized equivalent network, with same
                    % #nodes and degree distribution
                     connRbin=randomizer_bin_und(double(connAux2>0),1);
%                     % create GT metrics for subj s, cost c, randomisation r
                    globEffBinR(s,c,r)=efficiency_bin(connRbin);
                    locEffBinR(s,c,r)=mean(efficiency_bin(connRbin,1));
                    pathLengthBinR(s,c,r)=charpath(distance_bin(connRbin));
                    clustBinR(s,c,r)=mean(clustering_coef_bu(connRbin));
               end
        end
    end
end
           
% normalised pathlength (alpha), clustering coeff(gamma), small
% worldness(sigma)
% muGlobEffR=squeeze(mean(globEffR,3));
% muLocEffR=squeeze(mean(locEffR,3));
% muPathLR=squeeze(mean(pathLengthR,3));
% muClustR=squeeze(mean(clustR,3));
% alpha=pathLength./muPathLR;
% gamma=clust./muClustR;
% sigma=gamma./alpha;

muPathLBinR=squeeze(mean(pathLengthBinR,3));
muClustBinR=squeeze(mean(clustBinR,3));
alphaBin=pathLengthBin./muPathLBinR;
gammaBin=clustBin./muClustBinR;
sigmaBin=gammaBin./alphaBin; %smallworlness

save('GT_metricsSIFTadapted_costANDnormalized_130518')

%% Nodal degree group comparisons (90rois tested) 
% NodDegree = nrois x nsubs x vector network costs 
% You want to compare the groups for each node 
for s=1:Ntot
    for i=1:NumRois
%           LimcostintnodDeg(i,s)=mean(nodDegree(i,s,1:80));
%           LimcostintnodEffW(i,s)=mean(nodEffW(i,s,1:80));
%         LimcostintnodEffW_norm(i,s)=mean(nodEff_norm(i,s,1:80));
        LimcostintnodStrength_norm(i,s)=mean(nodStrength_norm(i,s,55:80));
       LimcostintnodStrength(i,s)=mean(nodStrength(i,s,55:80));
    end 
end
for i=1:NumRois
     PvaluesnodDeg(i)=anova1(LimcostintnodDeg(i,:), group, 'off');
     PvaluesnodEffW(i)=anova1(LimcostintnodEffW(i,:), group, 'off');
    PvaluesnodEffN(i)=anova1(LimcostintnodEffW_norm(i,:), group, 'off');
    PvaluesnodStrengthN(i)=anova1(LimcostintnodStrength_norm(i,:), group, 'off');
    PvaluesnodStrength(i)=anova1(LimcostintnodStrength(i,:), group, 'off');
end
[signif, pvalCrit]=fdr_bh(PvaluesnodDeg,0.05,'dep'); % nodal degrees
[signif2, pvalCrit2]=fdr_bh(PvaluesnodEffW,0.001,'dep'); % nodal trackweights (sum edges)
[signif3, pvalCrit3]=fdr_bh(PvaluesnodEffN,0.05,'dep'); % sum edges for normalized nodal efficiency
[signif4, pvalCrit4]=fdr_bh(PvaluesnodStrengthN,0.05,'dep');
[signif5, pvalCrit5]=fdr_bh(PvaluesnodStrength,0.01,'dep');
% factorPvals=0.05/pvalCrit2;correctedPval=Pvalues2*factorPvals;
% correctedPval(correctedPval>1)=1;

% Clustered barplot patients vs controls 
for i=1:NumRois
      [groupmeanNodEffs(i,:)]=[mean(LimcostintnodEffW(i,1:NsubsC)) mean(LimcostintnodEffW(i,NsubsC+1:NsubsC+NsubsD))];
      [groupmeanNodDeg(i,:)]=[mean(LimcostintnodDeg(i,1:NsubsC)) mean(LimcostintnodDeg(i,NsubsC+1:NsubsC+NsubsD))];
      [groupmeanNodEffN(i,:)]=[mean(LimcostintnodEffW_norm(i,1:NsubsC)) mean(LimcostintnodEffW_norm(i,NsubsC+1:NsubsC+NsubsD))];
      [groupmeanNodStrengthN(i,:)]=[mean(LimcostintnodStrength_norm(i,1:NsubsC)) mean(LimcostintnodStrength_norm(i,NsubsC+1:NsubsC+NsubsD))];
     [groupmeanNodStrength(i,:)]=[mean(LimcostintnodStrength(i,1:NsubsC)) mean(LimcostintnodStrength(i,NsubsC+1:NsubsC+NsubsD))];
     [groupStdNodStrength(i,:)]=[std(LimcostintnodStrength(i,1:NsubsC)) std(LimcostintnodStrength(i,NsubsC+1:NsubsC+NsubsD))];
end
figure; bh=bar(groupmeanNodStrength); 
nodThresh5=find(signif5==1);
title('Cost-integrated Nodal Strength by group', 'FontSize',12); 
ylabel('Nodal Strength','FontSize',20); 
xlabel('Node','FontSize',20); 
txt2 = [' * '] 
text(nodThresh5,repmat(max(max(groupmeanNodStrength)),size(nodThresh5,2),1),txt2, 'FontSize',20)
legend('Patients','Controls')

figure; bh=bar(groupmeanNodStrengthN); 
nodThresh4=find(signif4==1)
title('Cost-integrated Nodal Normalized Strength by group', 'FontSize',12); 
ylabel('Nodal Strength Normalized','FontSize',20); 
xlabel('Node','FontSize',20); 
txt2 = [' * '] 
text(nodThresh4,repmat(max(max(groupmeanNodEffN)),size(nodThresh4,2),1),txt2, 'FontSize',20)
legend('Patients','Controls')

% Order according to control nodal strength
[order1,order2]=sort(groupmeanNodStrength(:,2))
figure; bh=bar(groupmeanNodStrength(order2,:)); 
errorbar(groupmeanNodStrength(order2,:),groupStdNodStrength(order2,:)); 
for i=1:length(nodThresh5)
    nodThresh6(i)=find(order2==nodThresh5(i));
end
title('Cost-integrated Nodal Strength by group', 'FontSize',12); 
ylabel('Nodal Strength','FontSize',20); 
xlabel('Node','FontSize',20); 
txt2 = [' * '] 
text(nodThresh6,repmat(max(max(groupmeanNodStrength)),size(nodThresh6,2),1),txt2, 'FontSize',20)
legend('Patients','Controls')

%% Radiotherapy 
% Clustered barplots radiotherapy vs no radiotherapy vs controls 
RT=[0 1 0 1 1 1 1 0 0 0 0 1 0 1 1 1 1 1 1 1 1]
RT=RT(1,okC)
RTgroup=[RT, 2*ones(1,NsubsD)]
AaD=[3.2548 11.1068 7.2822 4.5836 10.7562 2.8603 18.4356 10.1808 ...
4.4822 7.7507 3.2466 13.2849 11.5616 8.7753 3.2438 7.1096 12.8685 ...
7.2137 4.8247 7.9562 13.9562]   
AaD=AaD(1,okC)

for i=1:NumRois
    [RTgroupmeanNodEff(i,:)]=[mean(LimcostintnodEffW(i,find(RTgroup==0))) mean(LimcostintnodEffW(i,find(RTgroup==1))) mean(LimcostintnodEffW(i,NsubsC+1:NsubsD))]
    [RTgroupmeanNodStrength(i,:)]=[mean(LimcostintnodStrength(i,find(RTgroup==0))) mean(LimcostintnodStrength(i,find(RTgroup==1))) mean(LimcostintnodStrength(i,NsubsC+1:NsubsD))]
    [RTgroupmeanNodStrengthN(i,:)]=[mean(LimcostintnodStrength_norm(i,find(RTgroup==0))) mean(LimcostintnodStrength_norm(i,find(RTgroup==1))) mean(LimcostintnodStrength_norm(i,NsubsC+1:NsubsD))]
end
%Within patient group comparison only 
for i=1:NumRois
    PvaluesRT(i)=anova1(LimcostintnodDeg(i,1:NsubsC), RT, 'off');
    Pvalues2RT(i)=anova1(LimcostintnodEffW(i,1:NsubsC), RT, 'off');
    PvaluesRTnodSt(i)=anova1(LimcostintnodStrength(i,1:NsubsC), RT, 'off');
    PvaluesRTnodStN(i)=anova1(LimcostintnodStrength_norm(i,1:NsubsC), RT, 'off');
end
[signif, pvalCrit]=fdr_bh(PvaluesRT,0.05,'dep');
nodThresh=find(signif==1);
[signif, pvalCrit]=fdr_bh(Pvalues2RT,0.05,'dep');
nodThresh=find(signif==1);
[signif, pvalCrit]=fdr_bh(PvaluesRTnodSt,0.05,'dep');
nodThresh=find(signif==1);
[signif, pvalCrit]=fdr_bh(PvaluesRTnodStN,0.05,'dep');
nodThresh=find(signif==1);

% Show mean for three subgroups 
[order1,order2]=sort(groupmeanNodStrength(:,2))
figure; bh=bar(RTgroupmeanNodStrength(order2,:)); 
legend('No Radiotherapy','Radiotherapy', 'Controls')
nodThresh5=find(Pvalues2RT<0.05)
title('Cost-integrated Nodal Strength by group', 'FontSize',12); 
ylabel('Nodal strength','FontSize',20); 
xlabel('Node','FontSize',20); 
txt2 = [' * '] 
text(nodThresh,repmat(max(max(RTgroupmeanNodStrength)),size(nodThresh6,2),1),txt2, 'FontSize',20)

% Show mean RT subgroups patients only 
for i=1:NumRois
    [RTgroupmeanPatNodStrength(i,:)]=[mean(LimcostintnodStrength(i,find(RTgroup==0))) mean(LimcostintnodStrength(i,find(RTgroup==1)))]
[RTgroupstdPatNodStrength(i,:)]=[std(LimcostintnodStrength(i,find(RTgroup==0))) std(LimcostintnodStrength(i,find(RTgroup==1)))]
end
figure; bh=bar(RTgroupmeanPatNodStrength(order2,:)); 
errorbar(RTgroupmeanPatNodStrength(order2,:),RTgroupstdPatNodStrength(order2,:)); 
title('Cost-integrated Nodal Strength by group', 'FontSize',12); 
ylabel('Nodal Strength','FontSize',20); 
xlabel('Node','FontSize',20); 
legend('No radiotherapy','Radiotherapy')

%% FA-based
nodStrengthFA=load('GT_FA','nodStrength')
nodStrengthFA=nodStrengthFA.nodStrength
nodEffWFA=load('GT_FA','nodEffW');nodEffWFA=nodEffWFA.nodEffW
for s=1:Ntot
    for i=1:NumRois
       LimcostintnodFAStrength(i,s)=mean(nodStrengthFA(i,s,55:80));
       LimcostintnodEffWFA(i,s)=mean(nodEffWFA(i,s,55:80));
    end 
end
for i=1:NumRois
      PvaluesnodFAStrength(i)=anova1(LimcostintnodFAStrength(i,:), group, 'off');
      PvaluesnodFAnodEffW(i)=anova1(LimcostintnodEffWFA(i,:), group, 'off');
end
[signif, pvalCrit]=fdr_bh(PvaluesnodFAnodEffW,0.05,'dep');
nodThresh=find(signif==1);
for i=1:NumRois
      [groupmeanNodFAStrength(i,:)]=[mean(LimcostintnodFAStrength(i,1:NsubsC)) mean(LimcostintnodFAStrength(i,NsubsC+1:NsubsD))];
 [groupmeanNodEffWFA(i,:)]=[mean(LimcostintnodEffWFA(i,1:NsubsC)) mean(LimcostintnodEffWFA(i,NsubsC+1:NsubsD))];
end
figure; bh=bar(groupmeanNodEffWFA); 
nodThresh=find(signif==1)
title('Cost-integrated Nodal Efficiency FA by group', 'FontSize',12); 
ylabel('Nodal Efficiency','FontSize',20); 
xlabel('Node','FontSize',20); 
txt2 = [' * '] 
text(nodThresh,repmat(max(max(groupmeanNodEffWFA)),size(nodThresh,2),1),txt2, 'FontSize',20)
legend('Patients','Controls')

%% NBS STATISTICS 
% Write tckweight mu-corrected connectomes for NBS 
Data_siftok=zeros(size(Data_all))
for i=1:size(Data_all,3)
    Data_siftok(:,:,i)=Data_all(:,:,i).*Data_all_mu(i)
end
Data_all=Data_siftok % Change your original tckweight connectome to corrected
indicator={'C01'; 'C02'; 'C04'; 'C05'; 'C06'; 'C07'; 'C08'; 'C09'; 'C10'; ...
    'C11'; 'C12'; 'C13'; 'C14'; 'C16'; 'C17'; 'C18'; 'C19'; 'C20'; 'C21'; ...
    'D01'; 'D02'; 'D03'; 'D04'; 'D05'; 'D06'; 'D07'; 'D08'; 'D09'; 'D10'; 'D11'; ...
    'D12'; 'D13'; 'D14'; 'D15'; 'D16'; 'D17'; 'D18'; 'D19'; 'D21'}
for i=1:length(indicator)
    name=[indicator{i} '_sift_connectome_AR_mu90Rois.csv'];
    dlmwrite(name,Data_siftok(1:90,1:90,i),'delimiter','\t');
end
% Write txtfile with list of connectomes listed for stats 
list={}; 
fileID = fopen('nbs_files.txt','w');
for i=1:length(indicator)
    list{i}=[indicator{i} '_sift_connectome_AR_mu90Rois.csv'];
    fprintf(fileID, '%s\n', list{i});
end 
fclose(fileID)
% Write design & contrast txt files 
list=[]; 
fileID = fopen('designMatrix.txt','w');
group2=[ones(1,NsubsC), zeros(1,NsubsD)]
for i=1:length(group)
    [list(i,:)]=[group(i) group2(i)];
    fprintf(fileID, '%s\n', num2str(list(i,:)));
end 
fclose(fileID); 

fileID = fopen('contrastPatLess.txt','w');
fprintf(fileID, '%s\n', '-1 1');
fclose(fileID);
