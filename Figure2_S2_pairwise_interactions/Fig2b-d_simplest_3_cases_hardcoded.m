% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 1/15/2021

% This script simulates the impact of one partitioning

%%
close all
clear all

addpath(genpath('../supporting functions'))
setFigDef

%% Basic simulation parameter setup
M = 2; % 10 populations
x0_orig = 1/M*ones(M,1); % uniform initial density across all populations
tend = 100; % the length of simulation

params_neutral{1}=[0.8;0.1];
params_neutral{2}=[0,0;0,0];
params_neutral{3}=[0,0;0,0];

params_neg{1}=[0.0;0.0];
params_neg{2}=[0,-1;0,0];
params_neg{3}=[0,0;0,0];

params_pos{1}=[1.5;0];
params_pos{2}=[0,0;0,0];
params_pos{3}=[0,+5;0,0];

params_pos2{1}=[0.8;0];
params_pos2{2}=[0,0;0,0];
params_pos2{3}=[0,+5;0,0];

%% set up the different segregation levels

cellTot = 1e2; % total number of cells
Ns = 6; % resolution of partitioning, number of partitionings to be simulated
part_max = cellTot*100; % the number of partitionings at the highest 
    % partitioning level
cell_desensity_max = 1e9; % the maximum number of cells in one population 
N_ = round(logspace(0,log10(part_max),Ns)); % the array of partitioning levels

vol_tot = 1; % total volume that is partitioned at each level

% sample the total number of cells into each local environment
BIsimpson_neu = zeros(10,Ns);
BIsimpson_neg = zeros(10,Ns);
BIsimpson_pos = zeros(10,Ns);
BIsimpson_pos2 = zeros(10,Ns);

for i = 1:10
    [nonEmpty, y0_] = seedInit(x0_orig,cellTot,N_,cell_desensity_max,vol_tot,1);

    figure(2)
    subplot(1,3,1)
    hold on
    yend_house = runSeg_2gamma(params_neutral, y0_, nonEmpty, cell_desensity_max, tend);
    [~,~, ~, BIsimpson_neu(i,:)] = plotBISeg(N_,yend_house,vol_tot*cell_desensity_max);
    plot(N_,BIsimpson_neu(i,:),'.','color',[0.7,0.7,0.7],'markersize',10)
    set(gca,'xscale','log')
    axis([1 part_max 1 2.1])

    subplot(1,3,2)
    hold on
    yend_house = runSeg_2gamma(params_neg, y0_, nonEmpty, cell_desensity_max, tend);
    [~,~, ~, BIsimpson_neg(i,:)] = plotBISeg(N_,yend_house,vol_tot*cell_desensity_max);
    plot(N_,BIsimpson_neg(i,:),'.','color',[0.7,0.7,0.7],'markersize',10)
    set(gca,'xscale','log')
    axis([1 part_max 1 2.1])
    set(gca,'xticklabel','')
    set(gca,'yticklabel','')    

    subplot(1,3,3)
    hold on
    yend_house = runSeg_2gamma(params_pos, y0_, nonEmpty, cell_desensity_max, tend);
    [~,~, ~, BIsimpson_pos(i,:)] = plotBISeg(N_,yend_house,vol_tot*cell_desensity_max);
    plot(N_,BIsimpson_pos(i,:),'.','color',[0.7,0.7,0.7],'markersize',10)
    set(gca,'xscale','log')
    axis([1 part_max 1 2.1])
    set(gca,'xticklabel','')
    set(gca,'yticklabel','')  
    
%     figure(2)
%     hold on
%     yend_house = runSeg_2gamma(params_pos2, y0_, nonEmpty, cell_desensity_max, tend);
%     [~,~, ~, BIsimpson_pos2(i,:)] = plotBISeg(N_,yend_house,vol_tot*cell_desensity_max);
%     plot(N_,BIsimpson_pos2(i,:),'.','color',[0.7,0.7,0.7],'markersize',10)
%     set(gca,'xscale','log')
%     axis([1 part_max 1 2.1])

end


%% plot mean and standard deviation
% Fig 2b-d
figure(2)
subplot(1,3,1)
hold on
plot(N_,mean(BIsimpson_neu),'-', 'LineWidth', 1,'color',[1,1,1]*0.4)
errorbar(N_,mean(BIsimpson_neu),std(BIsimpson_neu),'.', 'linewidth',1,'color',[1,1,1]*0.4)
plot(N_,mean(BIsimpson_neu),'ko', 'markersize',10, 'linewidth',2)
set(gca,'xscale','log')
xlabel('# of partitions')
ylabel('Simpson index')
axis square
% axis([1 part_max 1 2.1])


subplot(1,3,2)
hold on
plot(N_,mean(BIsimpson_neg),'k-', 'LineWidth', 1,'color',[1,1,1]*0.4)
errorbar(N_,mean(BIsimpson_neg),std(BIsimpson_neg),'k.', 'linewidth',1,'color',[1,1,1]*0.4)
plot(N_,mean(BIsimpson_neg),'ko', 'markersize',10, 'linewidth',2)
set(gca,'xscale','log')
axis square
% axis([1 part_max 1 2.1])


subplot(1,3,3)
hold on
plot(N_,mean(BIsimpson_pos),'k-', 'LineWidth', 1,'color',[1,1,1]*0.4)
errorbar(N_,mean(BIsimpson_pos),std(BIsimpson_pos),'k.', 'linewidth',1,'color',[1,1,1]*0.4)
plot(N_,mean(BIsimpson_pos),'ko', 'markersize',10, 'linewidth',2)
set(gca,'xscale','log')
axis square
% axis([1 part_max 1 2.1])

set(gcf,'position',[0 0 600 400])

% figure(2)
% plot(N_,mean(BIsimpson_pos2),'k-', 'LineWidth', 1,'color',[1,1,1]*0.4)
% errorbar(N_,mean(BIsimpson_pos2),std(BIsimpson_pos2),'k.', 'linewidth',1,'color',[1,1,1]*0.4)
% plot(N_,mean(BIsimpson_pos2),'ko', 'markersize',10, 'linewidth',2)
% set(gca,'xscale','log')
% 
% xlabel('# of partitions')
% ylabel('Simpson index')
% axis square
% set(gcf,'position',[0 0 200 400])
