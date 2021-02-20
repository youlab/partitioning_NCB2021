function [nonEmpty, Y0_]=seedInit(x0_orig,cellTot,N_,cell_desensity_max,vol_tot,FlagTakeOutZeros)

% x0_orig: a vector that dictates the relative abundance of each population
    % the sum of x0_orig should equal to 1

Ns = length(N_); % number of partitiong levels
M = length(x0_orig); % number of populations

Y0_ = cell(Ns,1); 
% initialize the output which is the sampled initial densities ready
% to be used as inputs for the simulations

nonEmpty = zeros(Ns,1); 
% initialize the vector that holds the number of non-empty local
% environments at each partitioning level


for i = 1:Ns % iterate through N partitioning levels
    tmpY0 = zeros(N_(i),M); 
    % hold the temporary sampled initial densities 
    % each row is one local community
    % each column is one population
    for j = 1:M % iterate through M populations

        tmpY0(:,j) = poissrnd(x0_orig(j).*cellTot/N_(i),N_(i),1)...
            /(cell_desensity_max*vol_tot/N_(i));
        % This samples the number of cells of one population in all
            % local environments belonging to the same partitioning level
        % x0_orig(j).*cellTot: the number of cells of this population in
            % total
        % x0_orig(j).*cellTot/numb_part(i): in each local
            % environment, the number of cells of this population on average
        % /(cell_desensity_max*vol_tot/numb_part(i)): 
            % scales the initial density based on the total carrying capacity
            % and the number of partitions
    end
    
if FlagTakeOutZeros==1
    Y0_{i} = tmpY0.*(tmpY0>(1/cell_desensity_max)); 

    % only keep the values that are greater than 1 cell

    inddd = find(sum(Y0_{i}(:,:),2)>=1e-9); 
        % find the local environments that have more than 1 cell
    Y0_{i} = Y0_{i}(inddd,:);
        % only keep the initial values of the local communities that are
        % non-empty 
    nonEmpty(i) = length(inddd);
        % count the number of local communities that are non-empty
else
    Y0_{i} = tmpY0;
end

end