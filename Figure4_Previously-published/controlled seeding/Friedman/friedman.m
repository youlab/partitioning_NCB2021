% 
% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 2/19/2021

% plotting data from
% Friedman, J., Higgins, L. M., & Gore, J. (2017). Community structure 
% follows simple assembly rules in microbial microcosms. 
% Nature ecology & evolution, 1(5), 1-7.

%% initialize environment

setFigDef

clear all

clc

%% data collected by manually counting
groupsize = [1,2,3,7,8];
presence = [8,8,8,5,3];

%%
figure(125)
subplot(1,3,1)
bar(presence,'facecolor',[1,1,1]*0.7)
xlabel('Local group size')
ylabel('Richness')
set(gca,'xtick',[1 2 3 4 5],'ytick',[0 2 4 6 8],'xticklabel',groupsize)
axis([0 6 0 8.5])
set(gca, 'XDir','reverse')
