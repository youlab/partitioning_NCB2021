% Created by Feilun Wu, feilunwu@gmail.com 
% & Yuanchi Ha, hayuanchi@gmail.com
% last edit: 08/29/2021

% This script generates remixing version of F2bcd
% the entire script can take more than 20 minutes (depending on the machine)

%% 
close all
clear all

%%
addpath(genpath('../supporting functions'))
setFigDef

%% 2c only negative + remixing process

M = 10; % 10 populations
x0_orig = 1/M*ones(M,1); % ones(M,1); %uniform initial density across all populations
tend = 100; % the length of simulation

cellTot = 1e3; % total number of cells
N_res = 6; % resolution of partitioning, number of partitionings to be simulated 6
part_max = cellTot*1; % the number of partitionings at the highest 10
    % partitioning level

frac_ = 1; % Only negative interactions
BIsimpsons_ = cell(1,1);
param_set = cell(5,1);
repeats = 10;
    
figure(5001)
BIsimpsons_ = zeros(repeats,N_res);
remixingTimes = 5;
remixingCount = 0;

for j = 1:repeats
    display(j)
    param_set{1} = 1; % [connectedness, neg_frac, max_delta, max_neg, max_pos]
    param_set{2} = frac_;
    param_set{3} = [0.0 0.0];
    param_set{4} = 0.8;
    param_set{5} = 0.0; 
    
    for r = 1:remixingTimes
        % The first run without remixing
        if remixingCount == 0
            %display(x0_orig)
            [N_, yend_house, ~, ~, nonEmpty] = getOneCurve(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e9,1);
            remixingCount = remixingCount + 1;
            yTemp = yend_house;
        else % Remix
            % Process each segregation level separately
            for k = 1: size(yTemp, 1)  
                % Sum each species across all environments from last run
                x0_orig = sum(yTemp{k}, 1)';
                % Normalize so that the sum of sums is 1
                x0_orig = normalize(x0_orig,'norm',1);   
                % Run as usual but only keep the result of corresponding
                % segregation level
                [N_, yend_house, ~, ~, nonEmpty] = getOneCurve(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e9,1);
                yTempFinal{k} = yend_house{k};
            end
            remixingCount = remixingCount + 1;
        end
    end
    yAllFinal{j} = yTempFinal;
    [~,~, ~, BIsimpsons_(j,:)] = plotBISeg(N_,yAllFinal{j},1e9);
    display('BIsimpsons: ')
    display(BIsimpsons_)
    %break
    hold on
    plot(N_,BIsimpsons_(j,:),'.','color',[0.7,0.7,0.7],'markersize',10)
    set(gca,'xscale','log')
    set(gca,'FontSize',20)
    axis([1 part_max 0.0 10])
end

plot(N_,mean(BIsimpsons_),'-','linewidth',1,'color',[1,1,1]*0.4)
plot(N_,mean(BIsimpsons_),'ko','markersize',10,'linewidth',2)

%% 2d only positive + remixing process

x0_orig = 1/M*ones(M,1);
frac_ = 0; % No neg interaction now
BIsimpsons_ = cell(1,1);
param_set = cell(5,1);
    
figure(5002)
BIsimpsons_ = zeros(repeats,N_res);
remixingTimes = 5;
remixingCount = 0;

for j = 1:repeats
    display(j)
    param_set{1} = 1; % [connectedness, neg_frac, max_delta, max_neg, max_pos]
    param_set{2} = frac_;
    param_set{3} = [0.0 2.0];
    param_set{4} = 0.0;
    param_set{5} = 3; 
    
    for r = 1:remixingTimes
        % The first run without remixing
        if remixingCount == 0
            %display(x0_orig)
            [N_, yend_house, ~, ~, nonEmpty] = getOneCurve(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e9,1);
            remixingCount = remixingCount + 1;
            yTemp = yend_house;
        else % Remix
            % Process each segregation level separately
            for k = 1: size(yTemp, 1)  
                % Sum each species across all environments from last run
                x0_orig = sum(yTemp{k}, 1)';
                % Normalize so that the sum of sums is 1
                x0_orig = normalize(x0_orig,'norm',1);   
                % Run as usual but only keep the result of corresponding
                % segregation level
                [N_, yend_house, ~, ~, nonEmpty] = getOneCurve(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e9,1);
                yTempFinal{k} = yend_house{k};
            end
            remixingCount = remixingCount + 1;
        end
    end
    yAllFinal{j} = yTempFinal;
    [~,~, ~, BIsimpsons_(j,:)] = plotBISeg(N_,yAllFinal{j},1e9);
    display('BIsimpsons: ')
    display(BIsimpsons_)
    %break
    hold on
    plot(N_,BIsimpsons_(j,:),'.','color',[0.7,0.7,0.7],'markersize',10)
    set(gca,'xscale','log')
    set(gca,'FontSize',20)
    axis([1 part_max 0.0 10])
end

plot(N_,mean(BIsimpsons_),'-','linewidth',1,'color',[1,1,1]*0.4)
plot(N_,mean(BIsimpsons_),'ko','markersize',10,'linewidth',2)


%% 2 both positive and negative + remixing process

x0_orig = 1/M*ones(M,1); % ones(M,1); %uniform initial density across all populations
frac_ = 0.5; % Half negative half positive interactions
BIsimpsons_ = cell(1,1);
param_set = cell(5,1);
    
figure(5003)
BIsimpsons_ = zeros(repeats,N_res);
remixingTimes = 5;
remixingCount = 0;

for j = 1:repeats
    display(j)
    param_set{1} = 1; % [connectedness, neg_frac, max_delta, max_neg, max_pos]
    param_set{2} = frac_;
    param_set{3} = [0.0 2.0];
    param_set{4} = 0.8;
    param_set{5} = 3; 
    
    for r = 1:remixingTimes
        % The first run without remixing
        if remixingCount == 0
            %display(x0_orig)
            [N_, yend_house, ~, ~, nonEmpty] = getOneCurve(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e9,1);
            remixingCount = remixingCount + 1;
            yTemp = yend_house;
        else % Remix
            % Process each segregation level separately
            for k = 1: size(yTemp, 1)  
                % Sum each species across all environments from last run
                x0_orig = sum(yTemp{k}, 1)';
                % Normalize so that the sum of sums is 1
                x0_orig = normalize(x0_orig,'norm',1);   
                % Run as usual but only keep the result of corresponding
                % segregation level
                [N_, yend_house, ~, ~, nonEmpty] = getOneCurve(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e9,1);
                yTempFinal{k} = yend_house{k};
            end
            remixingCount = remixingCount + 1;
        end
    end
    yAllFinal{j} = yTempFinal;
    [~,~, ~, BIsimpsons_(j,:)] = plotBISeg(N_,yAllFinal{j},1e9);
    display('BIsimpsons: ')
    display(BIsimpsons_)
    %break
    hold on
    plot(N_,BIsimpsons_(j,:),'.','color',[0.7,0.7,0.7],'markersize',10)
    set(gca,'xscale','log')
    set(gca,'FontSize',20)
    axis([1 part_max 0.0 10])
end

plot(N_,mean(BIsimpsons_),'-','linewidth',1,'color',[1,1,1]*0.4)
plot(N_,mean(BIsimpsons_),'ko','markersize',10,'linewidth',2)


