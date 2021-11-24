function [Y_end]=runSeg_2gamma_LV(params, Y0_, nonEmpty, carryingCap,tend)

% carryingCap

M = size(Y0_{1},2); % number of populations
options=odeset('NonNegative',1:M);
N_res = length(Y0_); % number of partitioning levels

Y_end = cell(N_res,1); % Initialize cell array

for i = 1:N_res
    Y_end{i} = zeros(nonEmpty(i),M); % each row is one local environment
    % each column is one population
end


%% run the simulation with various segregations

% older version
% yIndiv = (1-params{1}).*((1-params{1})>(1/cellTot));

% yIndiv = (1-params{1}); 
% yIndiv(yIndiv<=0)=0;

% there is a simple analytical solution when only one population is in a
% local environment
% segation levels

for i = 1:N_res
    for k = 1:nonEmpty(i)

%         if sum(Y0_{i}(k,:)>=(1/carryingCap))<=1 % only one or less strains are seeded in the local environment
%             Y_end{i}(k,:) = yIndiv'.*(Y0_{i}(k,:)>=(1/carryingCap));
%         else
%             
        [~,y] = ode45(@core_ode_lv,[0 tend],Y0_{i}(k,:), options, params{1}, params{2}, params{3});

        Y_end{i}(k,:) = y(end,:);
%         end
    end
%     i
end

fprintf('.')