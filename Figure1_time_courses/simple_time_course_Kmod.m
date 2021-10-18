% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 8/22/2021

% Simulating simple time courses
% Test the model construction, parameter generation, and the simulation process

%%
close all
clear all

addpath(genpath('../supporting functions'))
setFigDef


%% test pairwise communities

% Basic simulation parameter setup
N = 2; % 10 populations
x0 = 1/10*ones(N,1); % even initial density
tend = 100; % the length of simulation

% uncomment/comment to select the following pairwise communities  

% % no interactions
% connectedness = 1;
% neg_frac = 0.5;
% max_delta = 1;
% max_neg = 0;
% max_pos = 0;

% % only negative interaction
% connectedness = 1;
% neg_frac = 1;
% max_delta = 1;
% max_neg = 1;
% max_pos = 0;

% both negative and positive interaction
connectedness = 1;
neg_frac = 0.5;
maxmin_delta = [0 1.0];
max_neg = 1;
max_pos = 3;

% generate model parameters
params = param_generator(N,connectedness,neg_frac,...
    maxmin_delta,max_neg,max_pos);
params{4} = 0.8;

[t,y] = run_core_ode_Kmod(x0,tend,params);

figure(1)
plot(t,y)
xlabel('time (\tau)')
ylabel('density')


%%

% Basic simulation parameter setup
N = 10; % 10 populations
x0 = 1/1000*ones(N,1); % even initial density
tend = 200; % the length of simulation

% both negative and positive interaction
connectedness = 1;
neg_frac = 0.5;
maxmin_delta = [0 1.0];
max_neg = 1;
max_pos = 3;

% generate model parameters
for i = 1:5
params = param_generator(N,connectedness,neg_frac,...
    maxmin_delta,max_neg,max_pos,seeds(1));
params{4} = 1+(i-3)*0.1;

[t,y] = run_core_ode_Kmod(x0,tend,params);

figure(2)
subplot(1,5,i)
plot(t,y)
title(num2str(params{4}))
axis([0 200 0 1.2])
end

