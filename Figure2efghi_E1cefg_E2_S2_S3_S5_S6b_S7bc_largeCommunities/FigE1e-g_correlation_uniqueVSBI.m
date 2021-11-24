% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 1/15/2021

% This script generates extended data Figure 1e,f,g
% This script takes minutes to finish

%% 

close all
clear all

addpath(genpath('../supporting functions'))
setFigDef

%% Plot the number of unique types of communites

M = 10; % 2 populations
x0_orig = 1/M*ones(M,1); % uniform initial density across all populations
tend = 100; % the length of simulation

cellTot = 1e3; % total number of cells
N_res = 12; % resolution of partitioning, number of partitionings to be simulated 
part_max = cellTot*10; % the number of partitionings at the highest 
    % partitioning level
N_ = round(logspace(0,log10(part_max),N_res));


%% Extended data Figure 1f
cell_desensity_max = 1e9;
vol_tot = 1;
FlagTakeOutZeros = 1;

numberCount_all_house = zeros(100,N_res);
numberCount_eachCap_house = zeros(100,N_res);

for loopRand = 1:100 % generate 100 different sets of partitioning
    [nonempty, y0_]=seedInit(x0_orig,cellTot,N_,cell_desensity_max,vol_tot,FlagTakeOutZeros);

    numberCount_all = zeros(1,N_res);
    numberCount_eachCap = zeros(1,N_res);

    for i = 1:N_res
        pres = y0_{i}>(1/cell_desensity_max);
        numberCount_all(i) = size(unique(pres,'rows'),1);
        numberCount_eachCap(i) = mean(sum(pres,2));
    end
    numberCount_all_house(loopRand,:) = numberCount_all;

end
figure(343)
subplot(1,2,1)
plot(N_,numberCount_all_house,'.','color',[0.7,0.7,0.7],'markersize',10)
hold on
plot(N_,mean(numberCount_all_house),'-','color',[1,1,1]*0.4,'linewidth',1)
plot(N_,mean(numberCount_all_house),'ko','markersize',10)
set(gca,'yscale','linear','xscale','log')
xlabel('# of partitions')
ylabel({'Count of unique local',' communities at t_0'})
axis square


%% find for each population, how many unique communities it is sampled into
% Extended data Figure 1e, f

indiv_unique = zeros(N_res,M);

for i = 1:N_res
    pres = y0_{i}>(1/cell_desensity_max);

    for j = 1:M
        inds_tmp = y0_{i}(:,j)>(1/cell_desensity_max);
        pres = y0_{i}(inds_tmp,:)>(1/cell_desensity_max);
        indiv_unique(i,j) = size(unique(pres,'rows'),1);
    end

end

figure(343)

subplot(1,2,2)
plot(N_,indiv_unique,'.','color',[0.7,0.7,0.7],'markersize',10)
hold on
plot(N_,mean(indiv_unique,2),'-','color',[1,1,1]*0.4,'linewidth',1)
plot(N_,mean(indiv_unique,2),'ko','markersize',10)
set(gca,'yscale','linear','xscale','log')
axis square
xlabel('# of partitions')
ylabel({'Count containing' 'a population at t_0'})

set(gcf,'position',[0 0 600 300])

%% 
% plot the correlation between # of unique types of communities with 
% final biodiversity

% the parameters used to generate the .mat file
M = 10; % 2 populations
x0_orig = 1/M*ones(M,1); % uniform initial density across all populations
tend = 100; % the length of simulation
% % 
cellTot = 1e3; % total number of cells
N_res = 12; % resolution of partitioning, number of partitionings to be simulated
part_max = cellTot*10; % the number of partitionings at the highest 
    % partitioning level
N_ = round(logspace(0,log10(part_max),N_res));
% % 
cell_desensity_max = 1e9;
vol_tot = 1;
FlagTakeOutZeros = 1;


repeats = 10;

numSpec_ = [25,50,100];

BIsimpsons_ = cell(3,1);
numberCount_all_house_ = cell(3,1);

% figure(345)
for i = 1:3
    BIsimpsons_{i} = zeros(repeats,N_res);
    numberCount_all_house_{i} = zeros(repeats,N_res);
    M = numSpec_(i); % 2 populations
    x0_orig = 1/M*ones(M,1); 
    for j = 1:repeats
        % [connectedness, neg_frac, maxmin_delta, max_neg, max_pos]
        param_set = {1,0.5,[0 1.5],2,5}; 
        [N_, yend_house, y0_, ~] = getOneCurve(M,tend,param_set,x0_orig,...
            cellTot,N_res,part_max,cell_desensity_max,vol_tot);
        
        numberCount_all = zeros(1,N_res);
        numberCount_eachCap = zeros(1,N_res);

        for k = 1:N_res
            pres = y0_{k}>(1/cell_desensity_max);
            numberCount_all(k) = size(unique(pres,'rows'),1);
            numberCount_eachCap(k) = mean(sum(pres,2));
        end
        numberCount_all_house_{i}(j,:) = numberCount_all;    
        
        [~,~, ~, BIsimpsons_{i}(j,:)] = plotBISeg(N_,yend_house,cell_desensity_max*vol_tot);

    end

end

% % saving the output data
% save corr_unique_BI.mat numberCount_all_house_ BIsimpsons_
%%

figure(346)
for i = 1:3
    xthis_ = zeros(1,12*repeats);
    ythis_ = zeros(1,12*repeats);
    for j = 1:repeats
        hold on
        xthis_((12*(j-1)+1):12*j) = numberCount_all_house_{i}(j,:);
        ythis_((12*(j-1)+1):12*j) = BIsimpsons_{i}(j,:);
    end
    
    fittedParam = polyfit(log10(xthis_),ythis_,1);
    yPredict = polyval(fittedParam,log10(xthis_));
    plot(xthis_,yPredict,'color',[1,1,1]*0.28*i)    
    plot(xthis_,ythis_,'.','color',[1,1,1]*0.28*i,'markersize',10)
    
end
xlabel('Unique types of initial community')
ylabel('Simpson index')
set(gca,'xscale','log')
set(gca,'xtick',[1 10 100 1000])
axis('square')

