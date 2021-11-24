function [ts,ys] = sde_solve(y_init, tend, params)
%% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 7/18/2021

% stochastic simulator that can be used the same as 
% Simulating simple time courses with stochasticity using sde_solver.m
% Test the model construction, parameter generation, and the simulation
% process with sde_solver.m 

% y_init is a column vector
    %% Initialization and Utility
    dt = 0.1; % tested to ensure similarity to deterministic model
    sigma = 1e-2; % tested to ensure a meaningful noise level
%     sigma = 0; % for debugging...
    pd = makedist('Normal',0,sqrt(dt)); % Initialize the probability distribution for our 
                             % random variable with mean 0 and 
                             % stdev of sqrt(dt)

    dI = params{1};
    gamma = params{2};
    beta = params{3};

    ts    = 0 :dt :tend; % From t0-->t1 with N points
    ys    = zeros(length(y_init),length(ts));     % #species x # timepoints Matrix of zeros

    ys(:,1) = y_init;
    % find who are present in the initial seeding
    ys_present = y_init>1e-9;
    
    %% Computing the Process

        for i = 2:numel(ts)
            t = (i-1) .* dt;
            y = ys(:,i-1);
            
            dW = random(pd);

            ys(:,i) = y + (core_ode(t,y, dI,gamma,beta)).* dt + sigma .* dW;            

            ys(:,i) = ys(:,i).*(ys(:,i)>=0).* ys_present;
            % ensure ys is non-negative  
            
        end

        ys = transpose(ys);
        % transpose is to change the column vector, which is the output
        % of core_ode() to a row vector to store in ys. 
end

% S1: 
% make sure the euler method makes sense without noise; should be similar
% to the ode45 results

% S2: 
% introduce the noise in a seperate version 
% make sure the timestep is small enough that the results is converging
% then not change dt anymore 
% Only change the sigma term to scale dW

% Fig 2bcd
% Fig 2e 


