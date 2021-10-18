% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 1/15/2021

% This script generates Figure S2d
% the entire script can take more than 10 minutes (depending on the machine)

%% 

close all
clear all

addpath(genpath('../supporting functions'))
setFigDef


%% change proportions of negative interaction
% with a mixture of positive & negative interactions

M = 30; % 10 populations
x0_orig = 1/M*ones(M,1); % uniform initial density across all populations
tend = 100; % the length of simulation

cellTot = 1e3; % total number of cells
N_res = 6; % resolution of partitioning, number of partitionings to be simulated
part_max = cellTot*10; % the number of partitionings at the highest 
    % partitioning level

frac_ = [0.9, 0.7, 0.6, 0.5, 0.4, 0.1];
BIsimpsons_ = cell(6,1);
repeats = 10;

param_set = cell(5,1);
    
figure(20000)
for i = 1:6 % 6 fractions
    subplot(1,6,i)
    BIsimpsons_{i} = zeros(repeats,N_res);
    for j = 1:repeats
        param_set{1} = 1; % [connectedness, neg_frac, minmax_delta, max_neg, max_pos]
        param_set{2} = frac_(i);
        param_set{3} = [0.0 2]; %[0.2 4.2]; % [0 2]
        param_set{4} = 0.9; %0.7; % try 1
        param_set{5} = 0.9; %0.6; % try 0.6
        [N_, yend_house, ~, ~] = getOneCurve_LV(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e9,1);
        [~,~, ~, BIsimpsons_{i}(j,:)] = plotBISeg(N_,yend_house,1e9);
        display('BIsimpsons: ')
        display(BIsimpsons_)
        %break
        hold on
        plot(N_,BIsimpsons_{i}(j,:),'.','color',[0.7,0.7,0.7],'markersize',10)
        set(gca,'xscale','log')
        axis([1 part_max 0.0 30])
    end
    
    plot(N_,mean(BIsimpsons_{i}),'-','linewidth',1,'color',[1,1,1]*0.4)
    plot(N_,mean(BIsimpsons_{i}),'ko','markersize',10,'linewidth',2)
    %break
%     errorbar(N_,mean(BIsimpsons_{i}),std(BIsimpsons_{i}),'k.', 'LineWidth', 2)
end

figure(20000)
subplot(1,6,1)
xlabel('# of partitions')
ylabel('Simpson index')

for i = 1:6
    subplot(1,6,i)
    title((frac_(i)))
    if i>1
        set(gca,'xticklabel','')
        set(gca,'yticklabel','')   
    end
%     axis square
end

set(gcf,'position',[0 0 600 200])

%% change proportions of negative interaction
% with only positive interactions

M = 30; % 10 populations
x0_orig = 1/M*ones(M,1); % uniform initial density across all populations
tend = 100; % the length of simulation

cellTot = 1e3; % total number of cells
N_res = 6; % resolution of partitioning, number of partitionings to be simulated
part_max = cellTot*10; % the number of partitionings at the highest 
    % partitioning level

frac_ = [0.9, 0.7, 0.6, 0.5, 0.4, 0.1];
BIsimpsons_ = cell(6,1);
repeats = 10;

param_set = cell(5,1);
    
figure(20001)
for i = 1:6 % 6 fractions
    subplot(1,6,i)
    BIsimpsons_{i} = zeros(repeats,N_res);
    for j = 1:repeats
        param_set{1} = 1; % [connectedness, neg_frac, minmax_delta, max_neg, max_pos]
        param_set{2} = frac_(i);
        param_set{3} = [0.0 2]; 
        param_set{4} = 0.0;
        param_set{5} = 0.9;
        [N_, yend_house, ~, ~] = getOneCurve_LV(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e9,1);
        [~,~, ~, BIsimpsons_{i}(j,:)] = plotBISeg(N_,yend_house,1e9);
        display('BIsimpsons: ')
        display(BIsimpsons_)
        %break
        hold on
        plot(N_,BIsimpsons_{i}(j,:),'.','color',[0.7,0.7,0.7],'markersize',10)
        set(gca,'xscale','log')
        axis([1 part_max 0.0 30])
    end
    
    plot(N_,mean(BIsimpsons_{i}),'-','linewidth',1,'color',[1,1,1]*0.4)
    plot(N_,mean(BIsimpsons_{i}),'ko','markersize',10,'linewidth',2)
    %break
%     errorbar(N_,mean(BIsimpsons_{i}),std(BIsimpsons_{i}),'k.', 'LineWidth', 2)
end

figure(20001)
subplot(1,6,1)
xlabel('# of partitions')
ylabel('Simpson index')

for i = 1:6
    subplot(1,6,i)
    title((frac_(i)))
    if i>1
        set(gca,'xticklabel','')
        set(gca,'yticklabel','')   
    end
%     axis square
end

set(gcf,'position',[0 0 600 200])


%% change proportions of negative interaction
% with only negative interactions

M = 30; % 10 populations
x0_orig = 1/M*ones(M,1); % uniform initial density across all populations
tend = 100; % the length of simulation

cellTot = 1e3; % total number of cells
N_res = 6; % resolution of partitioning, number of partitionings to be simulated
part_max = cellTot*10; % the number of partitionings at the highest 
    % partitioning level

frac_ = [0.9, 0.7, 0.6, 0.5, 0.4, 0.1];
BIsimpsons_ = cell(6,1);
repeats = 10;

param_set = cell(5,1);
    
figure(20002)
for i = 1:6 % 6 fractions
    subplot(1,6,i)
    BIsimpsons_{i} = zeros(repeats,N_res);
    for j = 1:repeats
        param_set{1} = 1; % [connectedness, neg_frac, minmax_delta, max_neg, max_pos]
        param_set{2} = frac_(i);
        param_set{3} = [0.0 0]; 
        param_set{4} = 0.9;
        param_set{5} = 0.0;
        [N_, yend_house, ~, ~] = getOneCurve_LV(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e9,1);
        [~,~, ~, BIsimpsons_{i}(j,:)] = plotBISeg(N_,yend_house,1e9);
        display('BIsimpsons: ')
        display(BIsimpsons_)
        %break
        hold on
        plot(N_,BIsimpsons_{i}(j,:),'.','color',[0.7,0.7,0.7],'markersize',10)
        set(gca,'xscale','log')
        axis([1 part_max 0.0 30])
    end
    
    plot(N_,mean(BIsimpsons_{i}),'-','linewidth',1,'color',[1,1,1]*0.4)
    plot(N_,mean(BIsimpsons_{i}),'ko','markersize',10,'linewidth',2)
    %break
%     errorbar(N_,mean(BIsimpsons_{i}),std(BIsimpsons_{i}),'k.', 'LineWidth', 2)
end

figure(20002)
subplot(1,6,1)
xlabel('# of partitions')
ylabel('Simpson index')

for i = 1:6
    subplot(1,6,i)
    title((frac_(i)))
    if i>1
        set(gca,'xticklabel','')
        set(gca,'yticklabel','')   
    end
%     axis square
end

set(gcf,'position',[0 0 600 200])


