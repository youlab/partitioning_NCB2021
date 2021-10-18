% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 7/18/2021

% This script compares stochastic and deterministic simulations that 
% simulates the impact of multiple levels of partitioning on a
% community. 


%%
close all
clear all

addpath(genpath('../supporting functions'))
setFigDef

%%
% Basic simulation parameter setup
M = 10; % 10 populations
x0_orig = 1/M*ones(M,1); % uniform initial density across all populations
tend = 100; % the length of simulation

connectedness = 1;
neg_frac = 0.5;
minmax_delta = [0 1.5];
max_neg = 2;
max_pos = 5;

params = param_generator(M,connectedness,neg_frac,...
    minmax_delta,max_neg,max_pos,71);

%% set up the different segregation levels

cellTot = 1e3; % total number of cells
Ns = 5; % resolution of partitioning, number of partitionings to be simulated
part_max = cellTot*10; % the number of partitionings at the highest 
    % partitioning level
cell_desensity_max = 1e9; % the maximum number of cells in one population 
N_ = round(logspace(0,log10(part_max),Ns)); % the array of partitioning levels

vol_tot = 1; % total volume that is partitioned at each level

% sample the total number of cells into each local environment
[nonEmpty, y0_] = seedInit(x0_orig,cellTot,N_,cell_desensity_max,vol_tot,1);

%% run simulation to find the final density in each local environment

tic
yend_house_detm = runSeg_2gamma(params, y0_, nonEmpty, cell_desensity_max, tend);
toc

%% plot BI 
[bi,pres, densTot, BIsimpson] = plotBISeg(N_,yend_house_detm,cell_desensity_max*vol_tot);

figure(1)
hold on
% subplot(1,2,1)
plot(N_,BIsimpson,'o-')
set(gca,'xscale','log')


%%
tic
yend_house = runSeg_2gamma_Kmod(params, y0_, nonEmpty, cell_desensity_max, tend,0.5);
toc

%% plot BI
[bi,pres, densTot, BIsimpson] = plotBISeg(N_,yend_house,cell_desensity_max*vol_tot);

figure(1)
% subplot(1,2,2)
plot(N_,BIsimpson,'o-')
set(gca,'xscale','log')
legend('deterministic','randomize local carrying capacity')

%%
partLevel = 2;
figure(11)
subplot(1,2,1)
imagesc(yend_house{partLevel})
subplot(1,2,2)
imagesc(yend_house_detm{partLevel})
