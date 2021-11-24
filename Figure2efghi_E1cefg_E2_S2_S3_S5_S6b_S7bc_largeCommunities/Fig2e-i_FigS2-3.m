% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 1/15/2021

% This script generates Figure 2e-i, and Figure S2,3
% the entire script can take more than 10 minutes (depending on the machine)

%% 

close all
clear all

addpath(genpath('../supporting functions'))
setFigDef


%% Figure S2e
% change proportions of negative interaction

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
    
figure(13)
for i = 1:4
    subplot(5,4,i)
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

figure(13)
subplot(5,4,1)
%xlabel('# of partitions')
%ylabel('Simpson index')

for i = 1:4
    subplot(5,4,i)
    title((frac_(i)))
    if i>1
        set(gca,'xticklabel','')
        set(gca,'yticklabel','')   
    end
%     axis square
end

set(gcf,'position',[0 0 600 200])

%% Figure S2f
% Change number of species

numSpec_ = [5,10,20,50];
BIsimpsons_ = cell(4,1);
repeats = 10;
    
figure(13)
for i = 1:4
    subplot(5,4,i+4)
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

figure(13)
subplot(5,4,5)
%xlabel('# of species')
%ylabel('Simpson index')

for i = 1:4
    subplot(5,4,i+4)
    %title((frac_(i)))
    if i>1
        set(gca,'xticklabel','')
        %set(gca,'yticklabel','')   
    end
%     axis square
end

set(gcf,'position',[0 0 600 200])


%% Figure S2g
% Change the number of interactions


M = 20; % 20 populations
x0_orig = 1/M*ones(M,1); 
    
connectedness_ = [0.1,0.4,0.6,0.9];
BIsimpsons_ = cell(4,1);
repeats = 10;
    
figure(13)
for i = 1:4
    
    subplot(5,4,i+8)
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
figure(13)
subplot(5,4,9)
%xlabel('# of interactions')
for i = 1:4
    subplot(5,4,i+8)
    %title((frac_(i)))
    if i>1
        set(gca,'xticklabel','')
        set(gca,'yticklabel','')   
    end
%     axis square
end

%ylabel('Simpson index')

set(gcf,'position',[0 0 600 400])

%% Figure S2h
% change the strengths of interactions

M = 20; % 20 populations
x0_orig = 1/M*ones(M,1); 
    
strengths_ = [1, 1.5 , 2, 3];
BIsimpsons_ = cell(4,1);
repeats = 10;
    
figure(13)
for i = 1:4
    
    subplot(5,4,i+12)
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
subplot(5,4,13)
%xlabel('# of partitions')
%ylabel('Simpson index')
for i = 1:4
    subplot(5,4,i+12)
    %title((frac_(i)))
    if i>1
        set(gca,'xticklabel','')
        set(gca,'yticklabel','')   
    end
%     axis square
end

set(gcf,'position',[0 0 600 400])


%% Figure S2i
% change the ratio of the max strengths of interactions

M = 20; % 20 populations
x0_orig = 1/M*ones(M,1); 

strength = 2;
    
strengthSum = 1;
negs_ = [10/11, 0.5, 0.2, 1/11];
% [0.1, 0.4, 0.7, 0.9];

BIsimpsons_ = cell(4,1);
repeats = 10;
    
figure(13)
for i = 1:4
    subplot(5,4,i+16)
    BIsimpsons_{i} = zeros(repeats,N_res);
    
    for j = 1:repeats
        param_set{1} = 1; % [connectedness, neg_frac, minmax_delta, max_neg, max_pos]
        param_set{2} = 0.5;
        param_set{3} = [0 1.5];
        param_set{4} = negs_(i)*strength;
        param_set{5} = (strengthSum-negs_(i))*strength;
        display((strengthSum-negs_(i))/negs_(i))
        
        [N_, yend_house, ~, ~] = getOneCurve(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e9,1);
        [~,~, ~, BIsimpsons_{i}(j,:)] = plotBISeg(N_,yend_house,1e9);
        hold on
        plot(N_,BIsimpsons_{i}(j,:),'.','color',[0.7,0.7,0.7],'markersize',10)
        set(gca,'xscale','log')
        axis([1 part_max 0 20])
    end
    
    plot(N_,mean(BIsimpsons_{i}),'k-','linewidth',1,'color',[1,1,1]*0.4)
    plot(N_,mean(BIsimpsons_{i}),'ko','markersize',10,'linewidth',2)

    ratio = (strengthSum-negs_(i))/negs_(i);
    title(ratio)
    
end

figure(13)
subplot(5,4,17)
xlabel('# of partitions')
ylabel('Simpson index')
for i = 1:4
    subplot(5,4,i+16)
    %title((frac_(i)))
    if i>1
        set(gca,'xticklabel','')
        set(gca,'yticklabel','')   
    end
%     axis square
end

set(gcf,'position',[0 0 600 400])




%% Collect all Heatmap data -- can take up to 2 hours; load the heatmap directly instead
% Change number of species, interaction strengths, fraction of neg
% interactions

% Change number of species
numSpec_ = [5, 10, 20, 50, 100];
% Interaction strengths
strengths_ = [0.5, 1, 1.5, 2, 3, 4];
% Fraction of negative interactions
frac_ = [1.0, 0.9, 0.6, 0.4, 0.1, 0.0];

tend = 100; % the length of simulation

cellTot = 1e3; % total number of cells
N_res = 6; % resolution of partitioning, number of partitionings to be simulated
part_max = cellTot*10; % the number of partitionings at the highest 
    % partitioning level

BIsimpsons_Total = cell(4,4,4,1);
repeats = 10;

param_set = cell(5,1);
    
%figure(1000)
for s = 1:6 % strength of interactions
    param_set{4} = 0.4*strengths_(s);
    param_set{5} = 1.5*strengths_(s);
    
    for n = 1:5 % number of species
        M = numSpec_(n); % populations
        x0_orig = 1/M*ones(M,1); 
    
        for i = 1:6 % neg fractions
            %subplot(1,4,i)
            BIsimpsons_Total{s}{n}{i} = zeros(repeats,N_res);
            
            for j = 1:repeats
                param_set{1} = 1; % [connectedness, neg_frac, max_delta, max_neg, max_pos]
                param_set{2} = frac_(i);
                param_set{3} = [0 2];
                [N_, yend_house, ~, ~] = getOneCurve(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e9,1);
                [~,~, ~, BIsimpsons_Total{s}{n}{i}(j,:)] = plotBISeg(N_,yend_house,1e9);
                display('BIsimpsons: ')
                display(BIsimpsons_Total)
            end
        end
    end
end

%%
% Save the data here
save bigHeatmap.mat BIsimpsons_Total

%%
load bigHeatmap.mat

%% Create the heatmap for biodiversity
% Supplement Figure 3
fig = figure;
tiledlayout(6, 5, 'Padding', 'none', 'TileSpacing', 'compact'); 
pos = 1;
for s = 1:6
    for n = 1:5
        subplot(6,5,pos)
        pos = pos + 1;
        heatTable = [];
        x2 = [];
        for i = 1:6
            currentArray = BIsimpsons_Total{s}{n}{i};
            meanArray = mean(BIsimpsons_Total{s}{n}{i});
            meanArray = meanArray/numSpec_(n);
            heatTable = [heatTable; meanArray];
            % Fine max biodiversity for each role
            maximum = max(max(meanArray));
            [~, maxA] = find(meanArray==maximum);
            x2 = [x2 maxA];   
        end    
        
        % Plot the heatmap here
        colormap gray
        grid on;
        display(heatTable);
        h = imagesc(heatTable);
        %h.set(gca,'xgrid', 'on', 'ygrid', 'on')
        caxis manual
        %colorbar
        hold on

        % Plot the highlighter
%         y2 = [1 2 3 4 5 6];
%         plot(x2, y2, '-ro', 'MarkerSize', 8,...
%     'color', [1, 0.59, 0.55], 'MarkerFaceColor', [1, 0.59, 0.55]); %[0.53, 0.89, 0.31], [1, 0.59, 0.55]
        hold off

        if s == 6
            if n == 1
                xticks([1 6])
                yticks([1 2 3 4 5 6])
                xticklabels({'1', '1000'}) 
                yticklabels({'1.0', '0.9', '0.6', '0.4', '0.1', '0.0'})
                %[1.0, 0.9, 0.6, 0.4, 0.1, 0.0];
                xlabel('# of partitions') 
                ylabel('Fractions of negative interactions')                 
                cbh = colorbar;
                cbh.Position(1) = .95-cbh.Position(3);
                cbh.Position(2) = .9-cbh.Position(4)/2;
            else
                axis off;
            end
        else
            axis off;
        end
    end
    han=axes(fig,'visible','off'); 
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylh = ylabel(han,'Strength of interaction');
    ylh.Position(1) = ylh.Position(1) - 0.01;
    xlh = xlabel(han,'# of species');
    xlh.Position(2) = xlh.Position(2) - 0.01;
    xticklabels(han, {'5','10','20', '50', '100'})
end


%% Create the heatmap for biodiversity variance
fig = figure;
tiledlayout(6,5, 'Padding', 'none', 'TileSpacing', 'compact'); 
pos = 1;
for s = 1:6
    for n = 1:5
        subplot(6,5,pos)
        pos = pos + 1;
        heatTableV = [];
        x2 = [];
        for i = 1:6
            currentArray = BIsimpsons_Total{s}{n}{i};
            varArray = var(BIsimpsons_Total{s}{n}{i});
            varArray = varArray/numSpec_(n);
            heatTableV = [heatTableV; varArray];
        end
        
        colormap gray
        grid on;
        display(heatTableV);
        h = imagesc(heatTableV);
        %h.set(gca,'xgrid', 'on', 'ygrid', 'on')
        caxis manual
        %colorbar

        if s == 6
            if n == 1
                xticks([1 6])
                yticks([1 2 3 4 5 6])
                xticklabels({'1', '10000'}) 
                yticklabels({'1.0', '0.9', '0.6', '0.4', '0.1', '0.0'})
                xlabel('# of partitions') 
                ylabel('Fractions of negative interactions')                 
                cbh = colorbar;
                cbh.Position(1) = .95-cbh.Position(3);
                cbh.Position(2) = .9-cbh.Position(4)/2;
            else
                axis off;
            end
        else
            axis off;
        end
        %break

    end
    %break
    han=axes(fig,'visible','off'); 
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylh = ylabel(han,'Strength of interaction');
    ylh.Position(1) = ylh.Position(1) - 0.01;
    xlh = xlabel(han,'# of species');
    xlh.Position(2) = xlh.Position(2) - 0.01;
    xticklabels(han, {'5','10','20', '50', '100'})
end


%% Change hardcode in param generator for ideal results
% Supplement Figure 2

M = 2; % 2 populations
x0_orig = 1/M*ones(M,1); % ones(M,1); %uniform initial density across all populations
tend = 100; % the length of simulation

cellTot = 1e3; % total number of cells
N_res = 6; % resolution of partitioning, number of partitionings to be simulated 6
part_max = cellTot*1; % the number of partitionings at the highest 10


frac_ = 1; % Proportion of negative interactions
BIsimpsons_ = cell(1,1);
repeats = 10;

param_set = cell(5,1);
    
figure(100)

BIsimpsons_ = zeros(repeats,N_res);
for j = 1:repeats
    param_set{1} = 0; % [connectedness, neg_frac, max_delta, max_neg, max_pos]
    param_set{2} = frac_;
    param_set{3} = [0.0 2.0];
    param_set{4} = 0.0; 
    param_set{5} = 0; %3
    [N_, yend_house, ~, ~] = getOneCurve(M,tend,param_set,x0_orig,cellTot,N_res,part_max,1e9,1);
    display('size(yend_house): ')
    display((yend_house{1}))
    [~,~, ~, BIsimpsons_(j,:)] = plotBISeg(N_,yend_house,1e9);
    display('BIsimpsons: ')
    display(BIsimpsons_)
    %break
    hold on
    plot(N_,BIsimpsons_(j,:),'.','color',[0.7,0.7,0.7],'markersize',10)
    set(gca,'xscale','log')
    set(gca,'FontSize',16)
    axis([1 part_max 0.8 2.2])
end

plot(N_,mean(BIsimpsons_),'-','linewidth',1,'color',[1,1,1]*0.4)
plot(N_,mean(BIsimpsons_),'ko','markersize',10,'linewidth',2)
%legend('\delta_1 = 1.5','\delta_2 = 1.5')
    %break
%     errorbar(N_,mean(BIsimpsons_{i}),std(BIsimpsons_{i}),'k.', 'LineWidth', 2)


