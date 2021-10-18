% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 9/6/2021

% This script generates Extended Figure 2
% This script is used to investigate the distribution of biodiversity
% across local communities. 

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
N_res = 5; % resolution of partitioning, number of partitionings to be simulated
part_max = cellTot*10; % the number of partitionings at the highest 
    % partitioning level

frac_ = [1,0.5,0];
BIsimpsons_ = zeros(3,N_res);

param_set = cell(5,1);



for i = 1:3
    figure(21)
    subplot(1,3,i)
%     BIsimpsons_{i} = zeros(repeats,N_res);
        param_set{1} = 1; % [connectedness, neg_frac, minmax_delta, max_neg, max_pos]
        param_set{2} = frac_(i);
        param_set{3} = [0 2];
        param_set{4} = 0.8;
        param_set{5} = 3;
        rng(1)
        [N_, yend_house, ~, ~] = getOneCurve(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e9,1);
        [~,~, ~, BIsimpsons_(i,:)] = plotBISeg(N_,yend_house,1e9);
        hold on
        plot(N_,BIsimpsons_(i,:),'.','color',[0.7,0.7,0.7],'markersize',10)
        set(gca,'xscale','log')
        axis([1 part_max 0 20])
        if i == 1
            xlabel('# of partitions')
            ylabel('Simpson index')
        else
            set(gca,'xticklabel',[],'yticklabel',[])
%             yticklabel()
        end
        
    plot(N_,BIsimpsons_(i,:),'-','linewidth',1,'color',[1,1,1]*0.4)
    plot(N_,BIsimpsons_(i,:),'ko','markersize',10,'linewidth',2)
    
    localBI_ave = zeros(1,N_res);
    localBI_std = zeros(1,N_res);
    
    for j = 1:N_res
        BIs = simpsonInd_matrix(yend_house{j});
        count_zero = N_(j)-size(yend_house{j},1);
        BIs_addZeros = [BIs; zeros(count_zero,1)];
        
        localBI_ave(j) = mean(BIs_addZeros);
        localBI_std(j) = std(BIs_addZeros);
        
        [counts, bins] = hist(BIs_addZeros,10);
        figure(13+i)
        subplot(1,N_res,j)
        b = barh(bins,counts/N_(j));
        b.FaceColor = [1,1,1]*0.7;
        ylim([0 20])
        if j==1
            xlabel('normalized count')
            ylabel('Simpson index')
        else 
            
            set(gca,'xticklabel',[],'yticklabel',[])
%             yticklabel()
        end
    end
    
    figure(21)
    errorbar(N_,localBI_ave,localBI_std,'color',[0.5,0.5,1])
%     errorbar(N_,mean(BIsimpsons_{i}),std(BIsimpsons_{i}),'k.', 'LineWidth', 2)
end


