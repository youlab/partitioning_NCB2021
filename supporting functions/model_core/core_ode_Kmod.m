function dxdt = core_ode_Kmod(time,x,dI,gamma,beta,k)
% Adapted from the core model formulation
% with an addition of modulating carrying capacity (k)
% k is one value that directly represent the individual carrying capacity 
% of species. 

% gamma are negative values
% beta are positive values

% the ODE function
dxdt = x.*(1 - x./k + gamma*x)-(dI./(beta*x+1)).*x;

end
