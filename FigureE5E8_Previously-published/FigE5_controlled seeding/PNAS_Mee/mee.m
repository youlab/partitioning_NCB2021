% 
% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 2/19/2021

% plotting data from
% Mee, M. T., Collins, J. J., Church, G. M., & Wang, H. H. (2014). 
% Syntrophic exchange in synthetic microbial communities. 
% Proceedings of the National Academy of Sciences, 111(20), E2149-E2156.

%% initialize environment

setFigDef

clear all

clc
%% pairwise interactions (group size = 2)
[num,txt,~] = xlsread('pnas.1405641111.sd03.xlsx','2-member');


txt_strain_1 = txt(2:end,2); % strain 1
txt_strain_2 = txt(2:end,5);

fold_strain_1 = num(1:end,3); % fold change of strain 1
fold_strain_2 = num(1:end,6);

map_labels = {'M';'F';'K';'I';'R';'Y';'W';'T';'G';'C';'P';'L';'H';'S'};
map_labels_mat = cell2mat(map_labels);

%initialize--------------------------------------------------------

myMat2 = zeros(14); 
% myMat1 stores the fold change of each row strain, when it is paired with
% the strain indicated in the column.

for i = 1:91
    ind1 = strmatch(txt_strain_1{i},map_labels_mat);
    ind2 = strmatch(txt_strain_2{i},map_labels_mat);
    myMat2(ind1,ind2) = fold_strain_1(i);
    myMat2(ind2,ind1) = fold_strain_2(i);

end

relAbun2 = sum(myMat2,2)./sum(myMat2(:)); % sum all the pairs and calculate 
% relative abundance

%% 14 members (group size = 14)
[num14,txt,~]=xlsread('pnas.1405641111.sd03.xlsx','14&13-member mean');
relAbun14_1 = num14(:,4)./sum(num14(:,4)); % at day 3
relAbun14_2 = num14(:,5)./sum(num14(:,5)); % at day 4

%% 13 members (group size = 13)
% not plotted because not all 13 member communities were exausted


%% calculate biodiversity

relAbuns = [zeros(14,1) relAbun2 relAbun14_1]; 
% no strain can grow by themselves so the first column is populated with 0s

SimpsonBI = 1./sum(relAbuns.^2);

%% 3 group sizes as 3 levels of partitioning
groupsizes = [1,2,14];

%% plot subplot

figure(125)
subplot(1,3,2)
bar(SimpsonBI,'facecolor',[1,1,1]*0.7)
xlabel('Local group size')
ylabel('Simpson index')
set(gca,'xtick',[1 2 3 ],'xticklabel',groupsizes)
axis([0 4 0 8])
set(gca, 'XDir','reverse')
