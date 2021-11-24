function [N_, yend_house,y0_,params]=getOneCurve_Kmod(M,tend,param_set,x0_orig,cellTot,N_res,part_max,carrying_cap, vol_tot,Ksigma)
% Basic simulation parameter setup
% M: number of populations
% param_set: [connectedness, neg_frac, max_delta, max_neg, max_pos]

% x0_orig = 1/M*ones(M,1); % uniform initial density across all populations
% tend = 100; % the length of simulation

connectedness = param_set{1};
neg_frac = param_set{2};
minmax_delta = param_set{3};
max_neg = param_set{4};
max_pos = param_set{5};

params = param_generator(M,connectedness,neg_frac,minmax_delta,max_neg,max_pos);

%% set up the different segregation levels

% cellTot = 1e4; % total number of cells
% N_res = 5; % resolution of partitioning, number of partitionings to be simulated
% part_max = cellTot*10; % the number of partitionings at the highest 
%     % partitioning level
% cell_desensity_max = 1e9; % the maximum number of cells in one population 
N_ = round(logspace(0,log10(part_max),N_res)); % the array of partitioning levels

% vol_tot = 1; % total volume that is partitioned at each level

% sample the total number of cells into each local environment
[nonEmpty, y0_] = seedInit(x0_orig,cellTot,N_,carrying_cap,vol_tot,1);
% x0_orig,cellTot,N_,cell_desensity_max,vol_tot,FlagTakeOutZeros(1 is to
% take out 0s)

%% run simulation to find the final density in each local environment
tic
yend_house = runSeg_2gamma_Kmod(params, y0_, nonEmpty, carrying_cap, tend,Ksigma);
toc

% %% plot BI
% [bi,pres, densTot, BIsimpson] = plotBISeg(1:N_res,yend_house);
% 
% figure(11)
% hold on
% % subplot(2,2,1)
% plot(N_,BIsimpson,'o-')
% set(gca,'xscale','log')