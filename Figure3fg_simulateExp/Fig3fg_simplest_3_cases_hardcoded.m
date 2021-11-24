% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 9/5/2021

% This script simulates the experimental results of the 2 pairwise
% communities -- Fig3f,g

%%
close all
clear all

addpath(genpath('../supporting functions'))
setFigDef

%% Basic simulation parameter setup
M = 2; % 2 populations
x0_orig = 1/M*ones(M,1); % uniform initial density across all populations
tend = 100; % the length of simulation

params_neg{1}=[0;0];
params_neg{2}=[0,0;-0.7,0];
params_neg{3}=[0,0;0,0];

params_pos{1}=[0;1.1];
params_pos{2}=[0,-0.15;-0.15,0];
params_pos{3}=[0,0;0.8,0];


%% set up the different segregation levels

cellTot = 1e2; % total number of cells
Ns = 4; % resolution of partitioning, number of partitionings to be simulated
part_max = 10^3; % the number of partitionings at the highest 
    % partitioning level
cell_desensity_max = 1e9; % the maximum number of cells in one population 
N_ = [6 24 96 384]; % the array of partitioning levels

vol_tot = 1; % total volume that is partitioned at each level

% sample the total number of cells into each local environment
repeats = 10;

BIsimpson_neu = zeros(repeats,Ns);
BIsimpson_neg = zeros(repeats,Ns);
BIsimpson_pos = zeros(repeats,Ns);
BIsimpson_pos2 = zeros(repeats,Ns);



for i = 1:repeats
    [nonEmpty, y0_] = seedInit(x0_orig,cellTot,N_,cell_desensity_max,vol_tot,1);

    figure(1)


    subplot(1,2,1)
    hold on
    yend_house = runSeg_2gamma(params_neg, y0_, nonEmpty, cell_desensity_max, tend);
    [~,~, ~, BIsimpson_neg(i,:)] = plotBISeg(N_,yend_house,vol_tot*cell_desensity_max);
    plot(N_,BIsimpson_neg(i,:),'.','color',[0.7,0.7,0.7],'markersize',10)
    set(gca,'xscale','log')
    axis([2 part_max 1 2.1])
  

    subplot(1,2,2)
    hold on
    yend_house = runSeg_2gamma(params_pos, y0_, nonEmpty, cell_desensity_max, tend);
    [~,~, ~, BIsimpson_pos(i,:)] = plotBISeg(N_,yend_house,vol_tot*cell_desensity_max);
    plot(N_,BIsimpson_pos(i,:),'.','color',[0.7,0.7,0.7],'markersize',10)
    set(gca,'xscale','log')
    axis([2 part_max 1 2.1])
    
    


end


%% plot mean and standard deviation
figure(1)



subplot(1,2,1)
hold on
plot(N_,mean(BIsimpson_neg),'k-', 'LineWidth', 1,'color',[1,1,1]*0.4)
errorbar(N_,mean(BIsimpson_neg),std(BIsimpson_neg),'k.', 'linewidth',1,'color',[1,1,1]*0.4)
plot(N_,mean(BIsimpson_neg),'ko', 'markersize',10, 'linewidth',2)
set(gca,'xscale','log')
set(gca,'xtick',[10, 10^2, 10^3])
axis square
% axis([1 part_max 0 2.1])


subplot(1,2,2)
hold on
plot(N_,mean(BIsimpson_pos),'k-', 'LineWidth', 1,'color',[1,1,1]*0.4)
errorbar(N_,mean(BIsimpson_pos),std(BIsimpson_pos),'k.', 'linewidth',1,'color',[1,1,1]*0.4)
plot(N_,mean(BIsimpson_pos),'ko', 'markersize',10, 'linewidth',2)
set(gca,'xscale','log')
set(gca,'xticklabel','')
set(gca,'yticklabel','')  
axis square
% axis([1 part_max 0 2.1])

set(gcf,'position',[0 0 600 400])

