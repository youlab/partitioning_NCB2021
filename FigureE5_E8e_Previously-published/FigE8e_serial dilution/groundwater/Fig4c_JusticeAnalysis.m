%
% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 2/19/2021

% plotting data from
% Justice, N. B., Sczesnak, A., Hazen, T. C., & Arkin, A. P. (2017). 
% Environmental selection, dispersal, and organism interactions shape 
% community assembly in high-throughput enrichment culturing. 
% Applied and environmental microbiology, 83(20).

%% initialize environment

addpath(genpath('..\..\supporting functions'))
setFigDef

clear all
close all

clc

%% read in OTU sequencing data

otu_pool_raw = xlsread('otu_pool_full.xlsx');
otu_pool_raw = otu_pool_raw(:,2:end); % trim the first column, which is the indices

% trim the data table to leave only the OTUs with positive sequencing
% counts
ind_pos = sum(otu_pool_raw)>0;
OTUtot = sum(ind_pos); % OTUtot equals to 399, which is consistent with the number in
                        % the publication

otu_pool = otu_pool_raw(:,ind_pos);

% the number of OTU present in each pooled sample
numPres = sum(otu_pool>0,2);

%% normalize the OD to total of a pooled community = 1
tempMat = repmat(sum(otu_pool,2),1,OTUtot); % there are in total 1318 OTUs
otu_pool_norm = otu_pool./tempMat; % the 4th row of tempMat is 0 because there
                                    % was no sequences detected 


%% Find the relative abundance distribution of the original underground
% water sample.

binNumFit = 20;
relAbunBin = logspace(-7.5,-0.6,binNumFit); % define the two extremes of the 
% relative abundance

% =======================================================================
% This is the function of the histogram, ensuring a log-log linear
% relationship: 
% log(OTUcount) =  a*log(relAbunBin)+b 
syms a b
OTUcount = @(a,b) 10.^(a*log(relAbunBin)+b);

% The following define the two constrains for solving a and b
% 1st constraint: The total relative abundance should equal to 1.
% 2nd constraint: The total number of OTU is defined by "OTUtot".
% Note that OTU tot is greater than all the sequenced OTU, which is 399. 
OTUtot_orig = 5000;
constraints = [sum(OTUcount(a,b).*relAbunBin)==1, sum(OTUcount(a,b))==OTUtot_orig];
% =======================================================================

expStep = vpasolve(constraints,[a b]);
a_solved = expStep.a;
b_solved = expStep.b;
OTU_count_solved = round(OTUcount(a_solved,b_solved));

% check whether the sum is 1
% round(sum(OTU_count_solved.*relAbunBin))

% check the total number of OTU, this number should equal to OTUtot.
% sum(OTU_count_solved)

%% Plot the comparison between fitted distribution and experimental data

figure(1)
% plot fitted distribution with the final distribution
subplot(2,6,1)
plot(relAbunBin,OTU_count_solved,'.','color',[0.3,0.3,1]*0)
set(gca,'xscale','log','yscale','log')
legends = {'Estimate','10X','10^{2}X','10^{3}X','10^{4}X','10^{5}X'};
title(legends(1))
set(gca,'xscale','log','yscale','log','xticklabel','','yticklabel','')

hold on 
for i = 1:5
    subplot(2,6,i+1)

    sampleNum = i;
    hold on
    [a,b]=hist(log10(otu_pool_norm(sampleNum,otu_pool_norm(sampleNum,:)>0)),20);
    plot(10.^(b(a>0)),a(a>0),'.','color',[1,1,1]*0)
    set(gca,'xscale','log','yscale','log','xticklabel','','yticklabel','')
    title(legends(i+1))
    axis([1e-8 1e0 1 2000])
end


subplot(2,6,7)
plot(relAbunBin,OTU_count_solved,'.','color',[0.3,0.3,1]*0)
set(gca,'xscale','log','yscale','log')

hold on 
for i = 6:10
    subplot(2,6,i+2)
    sampleNum = i;
    hold on
    [a,b]=hist(log10(otu_pool_norm(sampleNum,otu_pool_norm(sampleNum,:)>0)),20);
    plot(10.^(b(a>0)),a(a>0),'.','color',[1,1,1]*0)
    set(gca,'xscale','log','yscale','log','xticklabel','','yticklabel','')
%     title(legends(i-4))
    axis([1e-8 1e0 1 2000])
end


subplot(2,6,7)
xlabel('Relative abundance')
ylabel('Count of OTU')

%% get the vector of an ideal community

cellDens_orig = 3.7e4; % cell density per ml of the original underground water sample
                        % according to the published paper
cellNumTot = cellDens_orig*1000; % assuming the first dilution samples were 
                            % sampled from 1000 ml of original sample,
                            % 3.55e7 cells. This number should be make sure
                            % the lowest relative abundance result in at
                            % least one cell. 


% initialize the index of OTU and the index of cells
indOTU = 0;
indCell = 0;

% initialize the output "indsCell_"
OTUtot_orig_actual = double(sum(OTU_count_solved));
indsOTU_ = 1:OTUtot_orig_actual;
indsCell_ = zeros(1,OTUtot_orig_actual); % stores the last index of the cells that belong
                                    % to the OTU. 
                                    
                                    
for i = 1:binNumFit % loop through all bins of relative abundances
    
    % define the indices of the first OTU and the last OTU that are at this
    % level of relative abundance. 
    OTU_first = indOTU+1;
    OTU_end = indOTU+OTU_count_solved(i);
    
    % number of cells for each OTU @ this level of relative abundance:
    num = round(relAbunBin(i)*cellNumTot);
    
    % numsTemp contains the final indices of cells that belong to 
    % each of all OTUs @ this level of relative abundance:
    numsTemp = (indCell+num):num:(indCell + OTU_count_solved(i)*num);
    
    % assign the cell indices to the output vector
    indsCell_(OTU_first:OTU_end) = numsTemp;

    % Update the end indices of OTU and Cells. 
    indOTU = OTU_end;
    indCell = indsCell_(OTU_end);

end


%% Sample from the original community
% Number of cells to be sampled: 
numCellSampled = round(3700*96*[1,0.1,0.01,0.001,0.0001]);

seeds_ =11:20;
repeats = 10; % repeat the sampling 10 times. 
sampled_ = cell(repeats,1);


% for k = 1:repeats
%     tic
%     sampled_tmp = zeros(5,OTUtot_orig_actual);
%     rng(seeds_(k))
%     for i = 1:5 % loop through all dilutions
%         
%         inds = randsample(cellNumTot_actual, numCellSampled(i)); % sample without replacement
%         inds_mat = repmat(inds,1,OTUtot_orig_actual);
%         cutoff_mat = repmat(indsCell_,numCellSampled(i),1);
% 
%         ID_mat = inds_mat-cutoff_mat; % this matrix helps identify which OTU a cell belongs to
%         for j = 1:numCellSampled(i)
%             OTUindTemp = find(ID_mat(j,:)<=0,1);
%             sampled_tmp(i,OTUindTemp) = sampled_tmp(i,OTUindTemp)+1;
%         end
%         
%     end
%     sampled_{k} = sampled_tmp;
%     toc
% end
% 
% save sampledOTU.mat sampled_

load sampledOTU.mat % uncomment the above loop to generate your own sampledOTU.mat

%%
richnessInit = zeros(5,repeats);
for i = 1:repeats
    richnessInit(:,i) = sum(sampled_{i}>=1,2);
end

richnessInit_ave = mean(richnessInit,2);
richnessInit_std = std(richnessInit,[],2);

%%

dilutions_ = logspace(1,5,5);

figure(14)
plot(dilutions_,numPres(1:5),'ro-')
hold on 
plot(dilutions_,numPres(6:10),'bo-')
% plot(richnessInit,'ko-')
errorbar(dilutions_,richnessInit_ave,richnessInit_std,'k','linewidth',2)

legend('NO_3','O_2','baseline')
set(gca,'yscale','log','xscale','log')
xlabel('Dilution')
ylabel('# OTU present')



%% Plot Figure 4c

figure(4234)
subplot(1,3,2)
hold on
% data points
plot(dilutions_,numPres(1:5)./richnessInit*100,'.','color',[1,1,1]*0.7,'markersize',12)
% error bar
errorbar(dilutions_,numPres(1:5)./richnessInit_ave*100,richnessInit_std./richnessInit_ave*100,'color',0.4*[1,1,1],'linewidth',1)
% connecting line
plot(dilutions_,numPres(1:5)./richnessInit_ave*100,'-','color',0.4*[1,1,1],'linewidth',1)
% marker -  average
plot(dilutions_,numPres(1:5)./richnessInit_ave*100,'ko','markersize',10)
axis square

set(gca,'xscale','log','xtick',[1e1 1e3 1e5])
xlabel('Dilution')
ylabel('% OTU present')
ytickformat(gca,'percentage')
axis([1 10^6 0 100])
title('NO_3')

subplot(1,3,3)
hold on
% data points
plot(dilutions_,numPres(6:10)./richnessInit*100,'.','color',[1,1,1]*0.7,'markersize',12)
% error bar
errorbar(dilutions_,numPres(6:10)./richnessInit_ave*100,richnessInit_std./richnessInit_ave*100,'color',0.4*[1,1,1],'linewidth',1)
% connecting line
plot(dilutions_,numPres(6:10)./richnessInit_ave*100,'-','color',0.4*[1,1,1],'linewidth',1)
% marker -  average
plot(dilutions_,numPres(6:10)./richnessInit_ave*100,'ko','markersize',10)
axis square

set(gca,'xscale','log','xticklabel','','yticklabel','')
axis([1 10^6 0 100])
title('O_2')

set(gcf,'position',[0 0 1000 500])

subplot(1,3,1)
plot(1,1,'.')
text(-0.5,0.5,'left empty')
text(-0.5,1,'intentionally')
axis off