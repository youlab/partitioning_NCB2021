function param = param_generator(N,connectedness,frac,minmax_del,max_neg,max_pos,rnd_seed)

% simplest function to generate the three sets of model parameters

% del, gam11, and gam22 all follow uniform distribution

% connectedness: the fraction of potential interactions that have
% interaction strengths (non-zero interaction strength)

% frac: the fraction of negative interaction over all interactions; This
% fraction is precise before zeroing out the diagnal of gam11 and gam22. 

% max_del: define the highest level of delta; minimum delta is 0.

% max_neg: define the maximum negative interaction strength; minimum
% negative interaction strength is 0

% max_pos: define the maximum positive interaction strength; minimum
% positive interaction strength is 0


% Set random number seed if a seed number is entered
if nargin==7
    rng(rnd_seed)
end

% Define the level of stress 
% the range is 0 to max_del, with an average at max_del/2
mindel = minmax_del(1);
maxdel = minmax_del(2);
del = rand(N,1)*(maxdel-mindel)+mindel; 

% Create a list of permutated indice for the interaction matrix
indsPerm = randperm(N*(N-1)); % exclude N diagonal values

% Initialize an interaction matrix
% An interaction is either positive or negative or neutral but cannot be
% both positive and negative. 
inter_init = zeros(N,N-1);

% Define the values of negative and positive interactions. 
neg_count = ceil(N*(N-1)*connectedness*frac);
pos_count = floor(N*(N-1)*connectedness*(1-frac));
gam = -rand(neg_count,1)*max_neg-0;
beta = +rand(pos_count,1)*max_pos+0;

% Assign the negative and positive interaction strengths to the interaction
% matrix, using the permutated indice. 
inter_init(indsPerm(1:neg_count)) = gam;
inter_init(indsPerm((neg_count+1):(pos_count+neg_count))) = beta;

% add the 0 diaganol values to the interaction matrix
inter = zeros(N,N);
for i = 1:N
    if i==1
        inter(i,2:N) = inter_init(i,1:N-1);
    elseif i==N
        inter(i,1:N-1) = inter_init(i,1:end);
    else
        inter(i,1:i-1) = inter_init(i,1:i-1);
        inter(i,i+1:N) = inter_init(i,i:end);
    end
end


% Split the negative interaction matrix and positive interaction matrix
gam_f = inter.*(inter<0); % negative interactions
beta_f = inter.*(inter>0); % positive interactions

% Wrap together the three sets of parameters
param = {del,gam_f,beta_f};

