function [numPerLevel] = seedInit_discrete_oneLevel(y0_orig,partitionLevel)

n = length(y0_orig); % number of species
cellTot = sum(y0_orig);
eachCell = repelem((1:n),y0_orig);

% sequential generation of the total number of cells in each partitioning
numPerLevel = zeros(1, partitionLevel);
cellLeft = cellTot;

for i = 1:(partitionLevel-1)
    lamb = cellLeft/(partitionLevel+1-i); % recompute lambda each time
    numPerLevel_temp = poissrnd(lamb);
    
    if numPerLevel_temp>cellLeft % sanity check
        numPerLevel_temp = cellLeft;
    end
    
    numPerLevel(i) = numPerLevel_temp;
    cellLeft = cellLeft-numPerLevel(i);
end

numPerLevel(end) = cellTot - sum(numPerLevel(1:(end-1)));

groupsIndx = repelem((1:partitionLevel),numPerLevel);

numPerLevel = createPartitionedVector(eachCell,groupsIndx,partitionLevel);

% y0_orig = 1./N;

% cellTot = 1E9;
% M = 9; % number of serial dilutions
% D = 5; % dilution factor
% D0 = 1E5;
% D_ = D0*D.^(1:M);% dilution

% M = length(D_); 
% Reps = D_/D0;
% 
% % RepsM = max(Reps);
% M = length(D_);
% N = length(y0_orig);
% % Y0_ = zeros(M,RepsM,N);
% Y0_ = cell(M,1);
% 
% % for i = 1:M
% %     Y0_{i}=zeros(Reps(i),N);
% % end
% 
% % Y0_(1,:,:) = repmat(y0_orig.*(y0_orig>1E-9),RepsM,1);
% 
% nonEmpty = zeros(M,1);
% 
% 
% for i = 1:M % different segregations
%     tmpY0 = zeros(Reps(i),N);
%     for j = 1:N
% %         size(poissrnd(y0_orig(j).*cellTot/D_(i),Reps(i),N)/cellTot*D_(i))
%         tmpY0(:,j) = poissrnd(y0_orig(j).*cellTot/D_(i),Reps(i),1)/cellTot*D_(i); % corrected for initial density
%     end
%     
% Y0_{i} = tmpY0.*(tmpY0>1E-9);
% inddd = find(sum(Y0_{i}(:,:),2)>=1e-4);
% Y0_{i} = Y0_{i}(inddd,:);
% nonEmpty(i) = length(inddd);

end