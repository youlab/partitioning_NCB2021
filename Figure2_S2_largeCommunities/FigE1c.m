% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 1/15/2021

% This script generates Extended Figure c

%% 

close all
clear all

addpath(genpath('../supporting functions'))
setFigDef

%%

M = 12; % 2 populations
x0_orig = 1/M*ones(M,1); % uniform initial density across all populations
tend = 100; % the length of simulation

cellTot = 5e2; % total number of cells
N_res = 8; % resolution of partitioning, number of partitionings to be simulated
part_max = 1e4; % the number of partitionings at the highest partitioning level

repeats = 10;
BIsimpsons_ = cell(2,1);

figure(12)

for i = 1:2
    
    BIsimpsons_{i} = zeros(repeats,N_res);
    if i==1
        param_set = {1,1,[0 0.2],1,0}; % [connectedness, neg_frac, max_delta, max_neg, max_pos]
        lightColor = [250,200,200]/255;
        solidColor = [250,18,18]/255;
    elseif i==2
        param_set = {1,0,[0 2],0,2};
        lightColor = [200,200,250]/255;
        solidColor = [25,25,250]/255;
    end
    
    for j = 1:repeats
        
        [N_, yend_house,~,~] = getOneCurve(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e7,1);
    %     M,tend,param_set,x0_orig,cellTot,N_res,part_max,carrying_cap, vol_tot
        [~,~, ~, BIsimpsons_{i}(j,:)] = plotBISeg(N_,yend_house,1e7);
        hold on
        p=plot(N_,BIsimpsons_{i}(j,:),'.','color',lightColor,'markersize',10);
        set(gca,'xscale','log')
        axis([1 part_max 2 12])
        
    end
    
    plot(N_,mean(BIsimpsons_{i}),'-','color',lightColor)
    errorbar(N_,mean(BIsimpsons_{i}),std(BIsimpsons_{i}),'LineWidth',1,'color',lightColor)
    plot(N_,mean(BIsimpsons_{i}),'o','color',solidColor,'markersize',10)
    
end

% save 12-member-X-figure.mat BIsimpsons_
%%
xlabel('# of partitioning')
ylabel('Simpson index')
axis('square')
