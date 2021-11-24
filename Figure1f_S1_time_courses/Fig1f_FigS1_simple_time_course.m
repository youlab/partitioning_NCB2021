% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 1/15/2021

% Simulating simple time courses
% Test the model construction, parameter generation, and the simulation process

%%
close all
clear all

addpath(genpath('../supporting functions'))
setFigDef


%% test pairwise communities

% Basic simulation parameter setup
N = 10; % 10 populations
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
max_neg = 0;
max_pos = 0.3;

% generate model parameters
params = param_generator_lv(N,connectedness,neg_frac,...
    maxmin_delta,max_neg,max_pos);

[t,y] = run_core_ode(x0,tend,params);

figure(1)
plot(t,y)
xlabel('time (\tau)')
ylabel('density')


%% Figure 
% Basic simulation parameter setup
N = 10; % 10 populations
x0 = 1/1000*ones(N,1); % even initial density
tend = 100; % the length of simulation

% both negative and positive interaction
connectedness = 1;
neg_frac = 0.5;
maxmin_delta = [0 2.0];
max_neg = 0.9; %1
max_pos = 0.9; %3

seeds = [36,51,65,91,68,99];

% generate model parameters
for i = 1:6
    params = param_generator(N,connectedness,neg_frac,...
        maxmin_delta,max_neg,max_pos,seeds(i));


    [t,y] = run_core_ode(x0,tend,params);

    figure(37)
    subplot(2,3,i)
    plot(t,y)
    %legend('y_1','y_2')
end

subplot(2,3,4)
xlabel('time (\tau)')
ylabel('density')


%% Figure 1e
% generate the main figure
% an example time course

N = 10;

connectedness = 1;
neg_frac = 0.5;
maxmin_delta = [0 1.5];
max_neg = 0.5;
max_pos = 5;

x0 = 1/10*ones(N,1); % even initial density

tend = 20;

kk=77; rng(kk)
params = param_generator(N,connectedness,neg_frac,...
    maxmin_delta,max_neg,max_pos);

[t,y] = run_core_ode(x0,tend,params);

paperColors = paperColor;
% colormap(paperColors)

figure(3)
set(gca, 'ColorOrder',paperColors, 'NextPlot','ReplaceChildren')
plot(t,y)
xlabel('Time')
ylabel('Density')


