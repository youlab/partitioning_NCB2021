%
% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 2/19/2021

% plotting data from
% Venturelli, O. S., Carr, A. V., Fisher, G., Hsu, R. H., Lau, R., Bowen, 
% B. P., ... & Arkin, A. P. (2018). 
% Deciphering microbial interactions in synthetic human gut microbiome communities. 
% Molecular systems biology, 14(6), e8157.


%% initialize environment

setFigDef

clear all

clc

%% monoculture OD, captured using GetDat Graph Digitizer

ER = 0.012;
EL = 0.115;
DP = 0.106;
FP = 0.118;
BH = 0.333;
CA = 0.382;
CH = 0.362;
PC = 0.182;
BO = 0.661;
BT = 0.737;
BU = 0.642;
BV = 0.745;
OD1 = [BH,CA,BU,PC,BO,BV,BT,EL,FP,CH,DP,ER];

groupsizes = [1,2,11,12];

%% read in data for pairwise interactions

[data2,txt2] = xlsread('inline-supplementary-material-2.xlsx');
% all 12 species
sp_names = {'BH','CA','BU','PC','BO','BV','BT','EL','FP','CH','DP','ER'};
% split names
pair_names = cell(66,2);

for i = 1:66
    this_pair_name = txt2{i,1};
    ind = find(isletter(this_pair_name),1);
    pair_names{i,1} = txt2{i,1}(ind:ind+1);
    pair_names{i,2} = txt2{i,1}(ind+2:ind+3);
end

myMat1 = zeros(12,12);

for i = 1:66
    ind1 = strmatch(pair_names{i,1},sp_names);
    ind2 = strmatch(pair_names{i,2},sp_names);
    myMat1(ind1,ind2) = data2(i,6);
    myMat1(ind2,ind1) = 1-data2(i,6);

end


%% read in data for 11-member communities and 12 member community
% 
[data11,txt11]=xlsread('inline-supplementary-material-4.xlsx');


%% calculate simpson index for each group size

simp_1 = 1/sum((OD1/sum(OD1)).^2);
simp_2 = 1/sum((sum(myMat1,2)/sum(sum(myMat1,2))).^2);
simp_11 = 1/sum((sum(data11(6:8:94,:))/12).^2);
simp_12 = 1/sum(data11(end-1,:).^2);

simps = [simp_1,simp_2,simp_11,simp_12];

%% plot Simpson's BI 
% for main figure 3

figure(125)
subplot(1,3,3)
bar(simps,'facecolor',[1,1,1]*0.7)
xlabel('Local group size')
ylabel('Simpson index')
set(gca,'xtick',[1 2 3 4],'xticklabel',groupsizes)
set(gca, 'XDir','reverse')

figure(125)
set(gcf,'position',[0 0 600 200])