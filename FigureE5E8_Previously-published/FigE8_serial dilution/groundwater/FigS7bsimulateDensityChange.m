% 
% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 2/19/2021

% generate data for FigS7b and plotting
% demonstrate that biodiversity stays fairly constant much when varying dilution
% alone


%% initialize the environment
addpath(genpath('..\..\..\supporting functions'))

setFigDef

clear all
close all

clc

%% Simulate count of present OTU with different starting dilutions

N = 1000;

connectedness = 1;
neg_frac = 0.5;
maxmin_delta = [0 1.5];
max_neg = 0.5;
max_pos = 2;

x0 = 1/1000*ones(N,1); % even initial density
dilutions_ = logspace(1,5,5);
tend = 200;

reps = 10;
richness = zeros(reps,5);
simpsonBI = zeros(reps,5);

% kk=77; rng(kk) 
for j = 1:reps
params = param_generator(N,connectedness,neg_frac,...
maxmin_delta,max_neg,max_pos);

    for i = 1:5
        [t,y] = run_core_ode(x0*dilutions_(i),tend,params);
        richness(j,i) = sum(y(end,:)>1e-7);
        simpsonBI(j,i) = simpsonInd(y(end,:));
    end

end

%%
figure(3)
subplot(1,2,1)
hold on
plot(dilutions_,richness,'.','color',[1,1,1]*0.7)
errorbar(dilutions_,mean(richness),std(richness),'color',[1,1,1]*0.4)
plot(dilutions_,mean(richness),'ko','markersize',12)

set(gca,'xscale','log','xtick',dilutions_)
xlabel('Dilution rate')
ylabel('Richness')
axis([10 1e5 0 200])

subplot(1,2,2)
hold on
plot(dilutions_,simpsonBI,'.','color',[1,1,1]*0.7)
errorbar(dilutions_,mean(simpsonBI),std(simpsonBI),'color',[1,1,1]*0.4)
plot(dilutions_,mean(simpsonBI),'ko','markersize',12)

set(gca,'xscale','log','xtick',dilutions_)
xlabel('Dilution rate')
ylabel('Simpson index')
axis([10 1e5 0 100])
