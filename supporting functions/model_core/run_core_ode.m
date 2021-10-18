function [t,y] = run_core_ode(x0,tend,params)

% this function solves the ODE model by specifiying the following values
% y0: the initial density
% tend: the end time point of the simulation
% params: the parameters of the model, including delta, 

N = length(x0);
options=odeset('NonNegative',1:N,'AbsTol',1e-9); 
% AbsTol is used to ignore x values lower than the density of 1 cell

[t,y] = ode45(@core_ode_lv,[0 tend],x0,options,params{1},params{2},params{3});

end