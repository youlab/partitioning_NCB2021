% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 1/15/2021

% This script generates figure S2a and S2b

%%
addpath(genpath('../supporting functions'))
setFigDef


clear all
close all

%% set up the parameters
M = 2;
names = {'1','2'};

% please run the script two times 
% one time with negative interaction and one time with positive interaction
% by commenting/uncommenting the following model parameter setup  

% ====== negative interaction ======
del = [0;0];
gam = [0,-1.5;0,0];
beta = [0,0;0,0];
indPlot = 1;
% ==================================

% ====== positive interaction ======
% del = [1.5;0];
% gam = [0,0;0,0];
% beta = [0,2;0,0];
% indPlot = 2;
% ==================================

gamma = gam+beta; % for plotting
params = {del, gam, beta};
x0_orig = 1/M*ones(1,M); 
N_ = [1,4,16];
cellTot = 20;
carryingCap = 1e5;
volumeTot = 1;

colorMapPaper = paperColor;
colorMapPaperLight = colorMapPaper([1:2:9],:);

tend = 200;
options = odeset('NonNegative',1:M);

%% sample initial densities and run simulations
rng(1) 
[nonEmpty,groupedInits] = seedInit(x0_orig,cellTot,N_,carryingCap,volumeTot,0);
yend_house = runSeg_2gamma(params, groupedInits, N_, carryingCap*volumeTot, tend);

%% No partitioning
figure(12345+indPlot)

k = 1;
subplot('Position',[0.05 0.5 0.3 0.3])
drawGraph_fixedLocation(gamma,names, 8*groupedInits{k}*1e3,1/(carryingCap/N_(k)), colorMapPaperLight,'circle'); 
addcircle(1.5)

subplot('Position',[0.05 0.1 0.3 0.3])
drawGraph_fixedLocation(gamma,names, yend_house{k}, 1/(carryingCap/N_(k)), colorMapPaperLight(:,:),'circle'); 
addcircle(1.5)

%% 2nd partitioning level, 4 partitions
k = 2;

rowNum = 1;
colNum = 1;

for i = 1:4 
    subplot('Position',[0.35 + 0.15*(rowNum-1) 0.5+0.17*(colNum-1) 0.15 0.15])
    drawGraph_fixedLocation(gamma, names, groupedInits{k}(i,:)*1e3, 1/(carryingCap/N_(k)), colorMapPaperLight(:,:),'circle'); 
    addcircle(1.5)

    subplot('Position',[0.35 + 0.15*(rowNum-1) 0.1+0.17*(colNum-1) 0.15 0.15])
    drawGraph_fixedLocation(gamma, names, yend_house{k}(i,:)/2, 1/(carryingCap/N_(k)), colorMapPaperLight(:,:),'circle'); 
    addcircle(1.5)
    
    colNum = colNum + 1;
    if colNum > 2
        colNum = 1;
        rowNum = rowNum+1;
    end
end

%% 3rd partitioning level, 16 partitions
k = 3;

rowNum = 1;
colNum = 1;

figure(12345+indPlot)
for i = 1:16 
    subplot('Position',[0.65 + 0.075*(rowNum-1) 0.5+0.085*(colNum-1) 0.075 0.075])
    drawGraph_fixedLocation(gamma, names,  groupedInits{k}(i,:)*1e3 ,1/(carryingCap/N_(k)), colorMapPaperLight(:,:),'circle'); 
    addcircle(1.5)

    subplot('Position',[0.65 + 0.075*(rowNum-1) 0.1+0.085*(colNum-1) 0.075 0.075])
    drawGraph_fixedLocation(gamma, names, yend_house{k}(i,:)/4,1/(carryingCap/N_(k)), colorMapPaperLight(:,:),'circle'); 
    addcircle(1.5)
    
    colNum = colNum + 1;
    if colNum > 4
        colNum = 1; 
        rowNum = rowNum+1;
    end
end

%% plot biodiversity index
[~,~,~,simpBI]=plotBISeg(N_,yend_house,carryingCap*volumeTot);
figure(1234)

subplot(1,2,indPlot)
plot(N_, simpBI,'ko','markersize',10,'linewidth',2)
hold on
plot(N_, simpBI,'-','linewidth',1,'color',0.4*[1,1,1])
set(gca,'xscale','log')
xlabel('Level of partitioning')
ylabel('Simpson index')
axis([1 16 1 2])

subplot(1,2,2)
xlabel('')
ylabel('')
set(gca,'xticklabel','','yticklabel','')
axis([1 16 1 2])