% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 1/15/2021

% This script makes sure the initial community spans over all 8 types of
% populations. 

% Extended Figure 9a-c


%% 

close all
clear all

addpath(genpath('../supporting functions'))
setFigDef


%% Setting up the basic range of simulation parameters

NN = 2; % resolution of the linear space of one param

resDel = 3;
resNeg = 3;
resPos = 2;

% % Default parameter ranges
linDel = repmat(linspace(0.5,1.5,NN),1,resDel); % stress the pops experience
gamNeg = repmat(-linspace(0,2,NN),1,resNeg); 
gamPos = repmat(linspace(0,5,NN),1,resPos);

% Total copies of the same archetype equals to: 
copies = resDel*resNeg*resPos;

[meshGamNeg,meshGamPos,meshDel] = meshgrid(gamNeg,gamPos,linDel);

% linearize each parameter
% the following values dictate the delta, the level of negative interaction
% and the level of positive interaction received by a population
del = meshDel(:);
linGamNeg = meshGamNeg(:); 
linGamPos = meshGamPos(:); 

M = length(del); 


%% generate interaction matrices 

% First create a matrix that combines the positive interaction parameters
% and the negative interaction parameters. 
% This combined matrix ensures that each strain is receiving the correct
% type of interaction from other species. 
gamNegHalf = ones(M,floor((M-1)/2)).*repmat(linGamNeg,1,floor((M-1)/2));
gamPosHalf = ones(M,ceil((M-1)/2)).*repmat(linGamPos,1,ceil((M-1)/2));
gam_noDiag = [gamNegHalf, gamPosHalf]; % this parameter contains all the values
% of the interaction matrix except the diagonal 

% Randomly permutate each row; Add 0 diagonal
% this permutation randomizes which population gives out a particular 
% interaction, without changing the type of interactions each population 
% receives 

rng(14)

gamma = zeros(M,M); 
for i = 1:(M-1)
    permedInteractionVector= gam_noDiag(i,randperm(size(gam_noDiag,2)));
    
    if i == 1
        gamma(i,:) = [0 permedInteractionVector];
    elseif i == M
        gamma(i,:) = [permedInteractionVector 0];
    else
        gamma(i,:) = [permedInteractionVector(1:i-1) 0 permedInteractionVector(i:end)];
    end
end


gammas = gamma.*(gamma<0);
betas = gamma.*(gamma>0);


%% set up the different partitioning levels

y0_orig = 6*ones(1,M); % 6 cells for each population
N_ = 6*4.^(0:1:5); % number of partitions


%% The initial density distribution

% an arbitrary number
originalCarryingK = 1e5; 

% cells per local environments
eachPartitionCarryingK = originalCarryingK./N_; 

% the threshold to determine whether a population is present at final time point
present_thresholds = 10./(eachPartitionCarryingK); 

options=odeset('NonNegative',1:M);
tend = 200;
params = {del, gamma.*(gamma<0), gamma.*(gamma>0)};


figure(45)
[pLocations,colorsEach]=drawGraph_customCircles(gamma, del, y0_orig, present_thresholds(1),0.07);

%% Run simulations

yEnd_House = zeros(length(N_),M);
yEnd_richness = zeros(1,length(N_));
yEnd_simpson = zeros(1,length(N_));

figure(323)
for i = 1:length(N_)
    
    % ================ initial conditions ========================
    groupedInit = seedInit_discrete_oneLevel(y0_orig,N_(i));
    y0 = groupedInit/eachPartitionCarryingK(i);
    % ============================================================
    
    yEnd = zeros(N_(i),M);

    for k = 1:N_(i)
        [t,y]=ode23(@core_ode,[0 tend],y0(k,:),options,params{1},params{2},params{3});
        yEnd(k,:) = y(end,:);    
    end
    
    subplot(1,length(N_),i)

    drawGraph_customCircles(gamma, del, sum(yEnd), present_thresholds(i),0.07,pLocations,colorsEach);

    yEnd_House(i,:) = sum(yEnd);
    yEnd_richness(i) = sum(sum(yEnd)>present_thresholds(1));
    yEnd_simpson(i) = simpsonInd(sum(yEnd));
    
end


%% quantify and plot the change of relative abundance of each archetype

typeNames = {'Fighter','Loser','Drama queen','Snowflake','Loner','Depressed','Exploiter','Follower'};

speciesParams =[linGamNeg linGamPos del];
[U,I] = unique(speciesParams,'row');
allIndices = zeros(length(del),8); % the grouping of the 8 types

% define the colors of each archetype based on parameters
color8types = cell(8,1);
for i = 1:8
    color8types{i,1} = [double(U(i,1)<mean(gamNeg)),double(U(i,2)>mean(gamPos)),double(U(i,3)<mean(linDel))];
    if sum(color8types{i,1})==3 % this is supposed to be white. change it to grey to visualize
        color8types{i,1}=[0.7,0.7,0.7];
    end
end

% total final density
totY_eachPartition = sum(yEnd_House,2); % total density of each partitioning level

% total pops grouped
groupedTot = zeros(length(N_),8); 

subplotInd = [6,5,8,7,2,1,4,3];

for i = 1:8 % looping through all 8 types
    figure(223)
    subplot(2,4,subplotInd(i))
    hold on
    lightC = (1-color8types{i})*0.6+color8types{i};

    [allIndices(:,i)]=ismember(speciesParams,U(i,:),'rows');
    for j = 1:length(N_)
        % adding the density of 18 copies of the same type
        groupedTot(j,i) = sum(yEnd_House(j,(allIndices(:,i)'==1)));
    end
    yEnd_norm_tmp = yEnd_House(:,(allIndices(:,i)'==1))./repmat(totY_eachPartition,1,copies);

    errorbar(N_,mean(yEnd_norm_tmp,2),std(yEnd_norm_tmp,[],2),'color',lightC,'LineWidth',1)
    plot(N_,yEnd_norm_tmp,'.','color',color8types{i});
    plot(N_,mean(yEnd_norm_tmp,2),'o','color',color8types{i},'LineWidth',2)
    
    if subplotInd(i)~=5
        set(gca,'xscale','log','xticklabels','','yticklabels','')
    else
        set(gca,'xscale','log')
    end
    axis([1 10000,0 0.04])
    title(typeNames{i})
    
    % plot Figure 5C
    figure(224)
    hold on
    errorbar(N_,mean(yEnd_norm_tmp,2),std(yEnd_norm_tmp,[],2),'color',lightC,'LineWidth',1)    
    plot(N_,mean(yEnd_norm_tmp,2),'o','color',color8types{i},'LineWidth',2,'markersize',10)
    set(gca,'xscale','log')
    
end

figure(223)
subplot(2,4,5)
xlabel('# of partitions')
ylabel('Relative abundance')

figure(224)
xlabel('# of partitions')
ylabel('Relative abundance')
axis([1 10000 0 0.03])


