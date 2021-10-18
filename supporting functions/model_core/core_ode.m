function dxdt = core_ode(time,x,dI,gamma,beta)
% the core model formulation

% gamma are negative values
% beta are positive values

% the ODE function
dxdt = x.*(1 - x + gamma*x)-(dI./(beta*x+1)).*x;

end
