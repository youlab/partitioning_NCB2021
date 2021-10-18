%% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 7/18/2021

% Simulating simple time courses with stochasticity using sde_solver.m
% Test the model construction, parameter generation, and the simulation
% process with sde_solver.m 

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


% change the following parameters to simulate different pairwise
% communities ===========
dI=[0.1;0.5];
gamma=[0,-1;0,0];
beta=[0,0;2,0];
% ===========

params = {dI,gamma,beta};

% run sde_solver
[ts,ys] = sde_solve(x0,tend,params);

% run ode45
[t,y] = run_core_ode(x0,tend,params);

figure(1)
subplot(1,2,1)
plot(ts,ys)
xlabel('time (\tau)')
ylabel('density')
title('stochastic')

subplot(1,2,2)
plot(t,y)
title('deterministic')
%% test large communities 

% Basic simulation parameter setup
N = 50; % # of populations
x0 = 1/10*ones(N,1); % even initial density
tend = 100; % the length of simulation

connectedness = 1;
neg_frac = 0.5;
maxmin_delta =[0,2];
max_neg = 2;
max_pos = 5;


% generate model parameters
params = param_generator(N,connectedness,neg_frac,...
    maxmin_delta,max_neg,max_pos);

% run sde_solver
[ts,ys] = sde_solve(x0,tend,params);

% run ode45
[t,y] = run_core_ode(x0,tend,params);

figure(2)
subplot(1,2,1)
plot(ts,ys)
xlabel('time (\tau)')
ylabel('density')
title('stochastic')

subplot(1,2,2)
plot(t,y)
title('deterministic')


%% regenerating the 6 timecourses presented in supplemental information using 
% stochastic simulations

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

seeds = [36,51,65,91,68,99];

% generate model parameters
for i = 1:6
params = param_generator(N,connectedness,neg_frac,...
    maxmin_delta,max_neg,max_pos,seeds(i));


[t,y] = run_core_ode(x0,tend,params);

figure(3)
subplot(2,3,i)
plot(t,y)

% run sde_solver
[ts,ys] = sde_solve(x0,tend,params);
figure(4)
subplot(2,3,i)
plot(ts,ys)

end

figure(3)
subplot(2,3,4)
xlabel('time (\tau)')
ylabel('density')
subplot(2,3,1)
title('deterministic')

figure(4)
subplot(2,3,4)
xlabel('time (\tau)')
ylabel('density')
subplot(2,3,1)
title('stochastic')

%% troubleshoot single timecourse stochastic vs deterministic comparisons
% 
% x0 = y0_{2}(3,:);
% tend = 100;
% 
% connectedness = 1;
% neg_frac = 0.5;
% minmax_delta = [0 1.5];
% max_neg = 2;
% max_pos = 5;
% 
% params = param_generator(M,connectedness,neg_frac,...
%     minmax_delta,max_neg,max_pos,71);
% 
% [t,y] = run_core_ode(x0,tend,params);
% 
% figure(123)
% subplot(1,2,1)
% plot(t,y)
% 
% % run sde_solver
% [ts,ys] = sde_solve(x0,tend,params);
% subplot(1,2,2)
% plot(ts,ys)