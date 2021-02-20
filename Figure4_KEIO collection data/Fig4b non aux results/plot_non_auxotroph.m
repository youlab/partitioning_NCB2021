% 
% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 2/19/2021

% process, analyze, and plot KEIO non-auxotroph community data

%% initialize environment

addpath(genpath('..\..\supporting functions'))
setFigDef

clear all
close all

clc

%% read in data

data = xlsread('non-auxotrophs.xlsx');
caliData = xlsread('non-auxotrophs.xlsx','even-mix calibration');

%% normalize abundance to calculate relative abundance
dataTemp = data./caliData(:,2);
dataTemp(isnan(dataTemp))=0;

dataTemp(dataTemp<0|dataTemp>1000)=0;

totReadCount = repmat(nansum(dataTemp),length(dataTemp),1);
realAbun = dataTemp./totReadCount;


%% calculate biodiversity index
BISimp = 1./nansum(realAbun.^2);


%% plot diversity against increasing partitioning levels
N_ = [6,24,96,384,1536]; % partitioning levels

figure(4234)
subplot(1,3,1)
hold on
plot(N_,BISimp(2:6),'ok','markersize',10)
plot(N_,BISimp(2:6),'-','markersize',10,'linewidth',1,'color',[1,1,1]*0.4)

set(gca,'xscale','log','xtick',[1e1 1e2 1e3])
axis([1,5000,35,50])
xlabel('# of partitions')
ylabel('Simpson index')

axis square

subplot(1,3,2)
plot(1,1,'.')
text (1,1,'left empty intentionally')
axis off
