% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 1/15/2021

% This script generate figures to demonstrate that spatial partitioning 
% implemented using controlled seeding also follows
% the same general principle

%% initialize environment

close all
clear all

addpath(genpath('../supporting functions'))
setFigDef



%% create random networks of different sizes, count how many have bell curve

M = 10; % Number of populations

kk_ = (61:70); % random number seeds

cellMax = 1e9;
D0 = 0.001;

options = odeset('NonNegative',1:M);
y0_orig = ones(1,M)/M;

tend = 100;

%% All competition
yend_house = cell(1,length(kk_));

connectedness = 1;
max_neg = 2;
max_pos = 1;

minmax_del_neg = [0 0]; 
frac_neg = 1; % 100% negative interactions

minmax_del_pos = [0 2]; 
frac_pos = 0; % 0% negative interactions (100% positive interactions)

bi = zeros(length(kk_),M,2);
pres = zeros(length(kk_),M,2);

lightColor1 = [250,200,200]/255;
lightColor2 = [200,200,250]/255;

solidColor1 = [250,18,18]/255;
solidColor2 = [25,25,250]/255;

figure(236)
for jjj = 1:length(kk_)
    
    yend_house{jjj} = cell(2,1);
    
    rng(kk_(jjj))
    [Y0_,nonEmpty] = seedInit_comb(y0_orig,cellMax,D0);
    
    params_neg = param_generator(M,connectedness,frac_neg,minmax_del_neg,max_neg,max_pos,kk_(jjj));
    yend_house{jjj}{1} = runSeg_2gamma(params_neg, Y0_, nonEmpty, cellMax,tend);
    [~,pres(jjj,:,1),~,bi(jjj,:,1)] = plotBISeg(1:M,yend_house{jjj}{1},cellMax*D0);
    
    params_pos = param_generator(M,connectedness,frac_pos,minmax_del_pos,max_neg,max_pos,kk_(jjj));
    yend_house{jjj}{2} = runSeg_2gamma(params_pos, Y0_, nonEmpty, cellMax,tend);
    [~,pres(jjj,:,2),~,bi(jjj,:,2)] = plotBISeg(1:M,yend_house{jjj}{2},cellMax*D0);
    
    subplot(1,4,1)
    hold on
    plot(pres(jjj,:,1),'.','color',lightColor1,'markersize',10)
    plot(pres(jjj,:,2),'.','color',lightColor2,'markersize',10)
    
    subplot(1,4,2)
    hold on
    plot(bi(jjj,:,1),'.','color',lightColor1,'markersize',10)
    plot(bi(jjj,:,2),'.','color',lightColor2,'markersize',10)

end

%%

figure(236)
subplot(1,4,1)
hold on
errorbar((1:M),mean(squeeze(pres(:,:,1))),std(squeeze(pres(:,:,1))),'LineWidth',1,'color',lightColor1)
errorbar((1:M),mean(squeeze(pres(:,:,2))),std(squeeze(pres(:,:,2))),'LineWidth',1,'color',lightColor2)
plot((1:M),mean(squeeze(pres(:,:,1))),'-','color',lightColor1,'linewidth',1)
plot((1:M),mean(squeeze(pres(:,:,2))),'-','color',lightColor2,'linewidth',1)
plot((1:M),mean(squeeze(pres(:,:,1))),'o','color',solidColor1,'markersize',6,'linewidth',2)
plot((1:M),mean(squeeze(pres(:,:,2))),'o','color',solidColor2,'markersize',6,'linewidth',2)

axis([0 11 0 12])
axis('square')
set(gca,'xtick',1:1:10)
ylabel('Richness')
xlabel('Level of partitioning')

subplot(1,4,2)
hold on
errorbar((1:M),mean(squeeze(bi(:,:,1))),std(squeeze(bi(:,:,1))),'LineWidth',1,'color',lightColor1)
errorbar((1:M),mean(squeeze(bi(:,:,2))),std(squeeze(bi(:,:,2))),'LineWidth',1,'color',lightColor2)
plot((1:M),mean(squeeze(bi(:,:,1))),'-','color',lightColor1,'linewidth',1)
plot((1:M),mean(squeeze(bi(:,:,2))),'-','color',lightColor2,'linewidth',1)
plot((1:M),mean(squeeze(bi(:,:,1))),'o','color',solidColor1,'markersize',6,'linewidth',2)
plot((1:M),mean(squeeze(bi(:,:,2))),'o','color',solidColor2,'markersize',6,'linewidth',2)


axis([0 11 0 12])
set(gca,'xtick',1:1:10,'xticklabel','')
ylabel('Simpson index')
axis('square')

%% both cooperation and competition
frac = 0.5; % 50% negative interactions
minmax_del = [0 2];
max_neg = 2;
max_pos = 5;


yend_house = cell(1,length(kk_));
bi = zeros(length(kk_),M);
pres = zeros(length(kk_),M);


figure(236)

for jjj = 1:length(kk_)
    
    rng(kk_(jjj))
    [Y0_,nonEmpty] = seedInit_comb(y0_orig,cellMax,D0);    
    params = param_generator(M,connectedness,frac,minmax_del,max_neg,max_pos,kk_(jjj));
    
    yend_house{jjj} = runSeg_2gamma(params, Y0_, nonEmpty, cellMax,tend);
    [~,pres(jjj,:),~,bi(jjj,:)] = plotBISeg(1:M,yend_house{jjj},cellMax*D0);
    
    subplot(1,4,3)
    hold on
    plot(pres(jjj,:),'.','color',[1,1,1]*0.7,'markersize',10)
  
    subplot(1,4,4)
    hold on
    plot(bi(jjj,:),'.','color',[1,1,1]*0.7,'markersize',10)

end

%%
subplot(1,4,3)
hold on
plot((1:M),mean(squeeze(pres(:,:,1))),'-','color',[1,1,1]*0.4,'LineWidth',1)
errorbar((1:M),mean(squeeze(pres(:,:,1))),std(squeeze(pres(:,:,1))),'LineWidth',1,'color',[1,1,1]*0.4)
plot((1:M),mean(squeeze(pres(:,:,1))),'o','color','k','markersize',6,'linewidth',2)

axis('square')
axis([0 11 0 12])
set(gca,'xtick',1:1:10,'xticklabel',1:1:10)
ylabel('Richness')
xlabel('Level of partitioning')

subplot(1,4,4)
hold on
plot((1:M),mean(squeeze(bi(:,:,1))),'-','color',[1,1,1]*0.4,'LineWidth',1)
errorbar((1:M),mean(squeeze(bi(:,:,1))),std(squeeze(bi(:,:,1))),'LineWidth',1,'color',[1,1,1]*0.4)
plot((1:M),mean(squeeze(bi(:,:,1))),'o','color','k','markersize',6,'linewidth',2)
axis('square')

axis([0 11 0 12])
set(gca,'xtick',1:1:10,'xticklabel','')
ylabel('Simpson index')
axis('square')