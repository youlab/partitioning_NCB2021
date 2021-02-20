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

M = 20; % 2 populations
x0_orig = 1/M*ones(M,1); % uniform initial density across all populations
tend = 100; % the length of simulation

cellTot = 1e3; % total number of cells
N_res = 6; % resolution of partitioning, number of partitionings to be simulated
part_max = cellTot*10; % the number of partitionings at the highest 
    % partitioning level

frac_ = [0.9,0.6,0.4,0.1];
BIsimpsons_ = cell(4,1);
repeats = 10;

param_set = cell(5,1);
    
figure(12)
for i = 1:4
    subplot(1,4,i)
    BIsimpsons_{i} = zeros(repeats,N_res);
    for j = 1:repeats
        param_set{1} = 1; % [connectedness, neg_frac, minmax_delta, max_neg, max_pos]
        param_set{2} = frac_(i);
        param_set{3} = [0 2];
        param_set{4} = 0.8;
        param_set{5} = 3;
        [N_, yend_house, ~, ~] = getOneCurve(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e9,1);
        [~,~, ~, BIsimpsons_{i}(j,:)] = plotBISeg(N_,yend_house,1e9);
        hold on
        plot(N_,BIsimpsons_{i}(j,:),'.','color',[0.7,0.7,0.7],'markersize',10)
        set(gca,'xscale','log')
        axis([1 part_max 0 20])
    end
    plot(N_,mean(BIsimpsons_{i}),'-','linewidth',1,'color',[1,1,1]*0.4)
    plot(N_,mean(BIsimpsons_{i}),'ko','markersize',10,'linewidth',2)
    
%     errorbar(N_,mean(BIsimpsons_{i}),std(BIsimpsons_{i}),'k.', 'LineWidth', 2)
end

figure(12)
subplot(1,4,1)
xlabel('# of partitions')
ylabel('Simpson index')

for i = 1:4
    subplot(1,4,i)
    title((frac_(i)))
    if i>1
        set(gca,'xticklabel','')
        set(gca,'yticklabel','')   
    end
%     axis square
end

set(gcf,'position',[0 0 600 200])

%% Figure S2d
% Change number of species

numSpec_ = [5,10,20,50];
BIsimpsons_ = cell(4,1);
repeats = 10;
    
figure(13)
for i = 1:4
    subplot(3,4,i)
    BIsimpsons_{i} = zeros(repeats,N_res);
    M = numSpec_(i); % 2 populations
    x0_orig = 1/M*ones(M,1); 
    for j = 1:repeats
        param_set{1} = 1; % [connectedness, neg_frac, minmax_delta, max_neg, max_pos]
        param_set{2} = 0.5;
        param_set{3} = [0 2];
        param_set{4} = 0.8;
        param_set{5} = 3;
        
        [N_, yend_house, ~, ~] = getOneCurve(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e9,1);
        [~,~, ~, BIsimpsons_{i}(j,:)] = plotBISeg(N_,yend_house,1e9);
        hold on
        plot(N_,BIsimpsons_{i}(j,:),'.','color',[0.7,0.7,0.7],'markersize',10)
        set(gca,'xscale','log')
        
    end
    plot(N_,mean(BIsimpsons_{i}),'k-','linewidth',1,'color',[1,1,1]*0.4)
    plot(N_,mean(BIsimpsons_{i}),'ko','markersize',10,'linewidth',2)

    title((numSpec_(i)))
end


%% Change the number of interactions


M = 20; % 20 populations
x0_orig = 1/M*ones(M,1); 
    
connectedness_ = [0.1,0.4,0.6,0.9];
BIsimpsons_ = cell(4,1);
repeats = 10;
    
figure(13)
for i = 1:4
    
    subplot(3,4,i+4)
    BIsimpsons_{i} = zeros(repeats,N_res);

    for j = 1:repeats
        param_set{1} = connectedness_(i); % [connectedness, neg_frac, minmax_delta, max_neg, max_pos]
        param_set{2} = 0.5;
        param_set{3} = [0 2];
        param_set{4} = 2;
        param_set{5} = 8;
        [N_, yend_house, ~, ~] = getOneCurve(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e9,1);
        [~,~, ~, BIsimpsons_{i}(j,:)] = plotBISeg(N_,yend_house,1e9);
        hold on
        plot(N_,BIsimpsons_{i}(j,:),'.','color',[0.7,0.7,0.7],'markersize',10)
        set(gca,'xscale','log')
        axis([1 part_max 0 20])
    end
    
    plot(N_,mean(BIsimpsons_{i}),'k-','linewidth',1,'color',[1,1,1]*0.4)
    plot(N_,mean(BIsimpsons_{i}),'ko','markersize',10,'linewidth',2)
    
    title((connectedness_(i)))
    
end


%% change the strengths of interactions

M = 20; % 20 populations
x0_orig = 1/M*ones(M,1); 
    
strengths_ = [1,1.5,2,3];
BIsimpsons_ = cell(4,1);
repeats = 10;
    
figure(13)
for i = 1:4
    
    subplot(3,4,i+8)
    BIsimpsons_{i} = zeros(repeats,N_res);

    for j = 1:repeats
        param_set{1} = 1; % [connectedness, neg_frac, minmax_delta, max_neg, max_pos]
        param_set{2} = 0.5;
        param_set{3} = [0 1.5];
        param_set{4} = 0.4*strengths_(i);
        param_set{5} = 1.5*strengths_(i);
        
        [N_, yend_house, ~, ~] = getOneCurve(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e9,1);
        [~,~, ~, BIsimpsons_{i}(j,:)] = plotBISeg(N_,yend_house,1e9);
        hold on
        plot(N_,BIsimpsons_{i}(j,:),'.','color',[0.7,0.7,0.7],'markersize',10)
        set(gca,'xscale','log')
        axis([1 part_max 0 20])
    end
    plot(N_,mean(BIsimpsons_{i}),'k-','linewidth',1,'color',[1,1,1]*0.4)
    plot(N_,mean(BIsimpsons_{i}),'ko','markersize',10,'linewidth',2)

    title((strengths_(i)))
    
end

figure(13)
subplot(3,4,9)
xlabel('# of partitions')
ylabel('Simpson index')

set(gcf,'position',[0 0 600 400])
