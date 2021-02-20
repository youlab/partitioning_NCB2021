function [Y0_, typeCount]=seedInit_comb(y0_orig,cellTot,D0)

% y0_orig = 1./N; % this is the fraction of each type of cells

% D0 is the inital total cell density

M = length(y0_orig); % M is the number of populations 
allComb = dec2bin(2^M-1:-1:0)-'0'; % creates all combinations of M strains
% "1" indicates the strain is present
% "0" indicates the strain is not present

gs = sum(allComb,2); % calculate the group size of each combination

Y0_ = cell(M,1); % initialize the output, which are the initial conditions 


typeCount = zeros(M,1); % count the types of initial community at each 
    % group size


for i = 1:M
    % take out all the combinations that have group size equal to "i"
    % "gs" means group size
    ind_gs = allComb(gs==i,:); 
    
    % calculate the initial density
    % Density of each strain is normalized by "i" to ensure the total
    % density equals to D0 regardless of gs value
    Y0_{M+1-i,1} = ind_gs/i*D0; 
    
    typeCount (M+1-i) = size(Y0_{M+1-i,1},1);
end


end