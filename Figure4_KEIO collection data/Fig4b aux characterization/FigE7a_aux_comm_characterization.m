% 
% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 2/19/2021

% Plot OD values of the auxotrophic community at 48 hrs 
% 6 different [CA] were used to culture the community
% 7 technical replicates
% Extended Figure 7a

%% initialize environment

addpath(genpath('..\..\supporting functions'))
setFigDef

clear all
close all

clc


%% read in data

data = xlsread('48 aux char -2-16-20.xlsx', 'Summary');

%%

meanOD = data(2:7,:);
stdvOD = data(11:16,:);

dilutions = data(1,:);
initialDensity = 1.4e9./dilutions/5;

%%
figure(1)
hold on
for i= 1:6
    set(gca,'xscale','log')
    plot(initialDensity,meanOD(i,:),'o','color',[1,1,1]*0.12*i,'markersize',8)
    errorbar(initialDensity',meanOD(i,:)',stdvOD(i,:)','.','color',[1,1,1]*0.15*i)
    plot(initialDensity,meanOD(i,:),'-','linewidth',1,'color',[1,1,1]*0.15*i)
    
end

ylabel('OD')
xlabel('initial density (CFU/well)')