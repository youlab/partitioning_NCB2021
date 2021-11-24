% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 1/15/2021

% This script plot the data from the sponge experiment
% Fig 5d

%% initialize environment

close all
clear all

addpath(genpath('..\supporting functions'))
setFigDef

%% read in data
data = xlsread('sponge_results.xlsx'); 

%% normalize by even mix and calculate relative abundance
evenMix = data(:,7);
spongeExp = data(:,1:6);

spongeExp_norm = spongeExp./(evenMix);
spongeExp_norm_RA = spongeExp_norm./nansum(spongeExp_norm);

%% calcualte simpson index and plot

sampleCount = 6;
simp_sponge = zeros(sampleCount,1);

for i = 1:sampleCount
    simp_sponge(i) = simpsonInd(spongeExp_norm_RA(:,i));
end

simp_sponge_mat_tmp = reshape(simp_sponge,3,2);
simp_sponge_mat = simp_sponge_mat_tmp([1,3,2],:);

figure(1)
hb = bar(simp_sponge_mat);
hb(1).FaceColor = [1,1,1]*0.8;
hb(2).FaceColor = [1,1,1]*0.4;
xticklabels({'0.001%', '0.005%', '0.02%'})
ylabel('Simpson index')
pos = [0 0 280 200];
set(gcf,'position',pos)
legend({'-','+'})

