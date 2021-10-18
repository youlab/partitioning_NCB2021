function dxdt = core_ode(time, x, dI, gamma, beta) %dI,
% the core model formulation

% gamma are negative values
% beta are positive values

%%% gLV version with intrinsic cap -- supplement formulation
dxdt = x.*(1 - x - dI + gamma*x + beta*x).*(1-x);

%%% gLV version with a hard cap
% x(x>1)=1;
% dxdt = x.*(1 - x - dI + gamma*x + beta*x); % LV
% setIndex = (x>=1);
% finalIndex = zeros(size(dxdt));
% 
% for i = 1:length(dxdt)
%     if setIndex(i) == 1
%         if dxdt(i) > 0
%             finalIndex(i) = 1;
%         end
%     else
%         finalIndex(i) = 0;
%     end
%     
% end
% 
% dxdt(finalIndex == 1) = 0;


end
