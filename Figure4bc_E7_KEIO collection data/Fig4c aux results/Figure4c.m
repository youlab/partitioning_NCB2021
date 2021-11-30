% 
% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 2/19/2021

% process, analyze, and plot KEIO auxotroph community data
% Fig 4c

addpath(genpath('..\..\supporting functions'))
setFigDef

clear all
close all

clc

%% read in data for 1st technical replicate and normalize 

% read in data
data_1 = xlsread('auxotroph_data.xlsx','Rep1'); 

% normalize by even mix
controlInd = 92; % the index of even mix control sample
data_1_norm = data_1(:,:)./data_1(:,controlInd);
data_1_norm(isnan(data_1_norm))=0;
data_1_norm((data_1_norm==Inf))=0;
data_1_norm_RA = data_1_norm./sum(data_1_norm); % normalized relative abundance


%% read in data for 2nd technical replicate and normalize

% read in data
data_2 = xlsread('auxotroph_data.xlsx','Rep2'); 

% normalize by even mix
data_2_norm = data_2(:,:)./data_2(:,controlInd);
data_2_norm(isnan(data_2_norm))=0;
data_2_norm((data_2_norm==Inf))=0;
data_2_norm_RA = data_2_norm./sum(data_2_norm); % normalized relative abundance


%% calculate Simpson's Index 

sampleCount = 92;

simp_1_pre = zeros(1,sampleCount);
simp_2_pre = zeros(1,sampleCount);

for i = 1:sampleCount
    
    simp_1_pre(i) = simpsonInd(data_1_norm_RA(:,i));
    simp_2_pre(i) = simpsonInd(data_2_norm_RA(:,i));

end

%% Reshape and rearrange the data for easier downstream processing
% add padding for easy reshaping
simp_1_pre = [simp_1_pre, 0,0,0,0];
simp_2_pre = [simp_2_pre, 0,0,0,0];

simp_1 = sampleRearrange(simp_1_pre);
simp_2 = sampleRearrange(simp_2_pre);

%% ================ Plot main figure ====================================
% and calculate statistical significance of trends

% 0% [CA] is not plotted because the high variability due to overall low cell 
% count

xaxisValues = [6,24,96,384,1536];
modelX = cell(5,1); % initial a cell to store  fitted models for the five 
% [CA] concentrations

titles = {'0.0002%','0.001%','0.005%','0.02%','0.1%'}; % [CA]

slopeX = zeros(7,2); % housekeeping variable
% first column is slope, second column is P valued

% 1st row: [CA] = 0.0002%
% 2nd row: [CA] = 0.001%
% 3rd~5th rows: [CA] = 0.005%
% 6th row: [CA] = 0.02%
% 7th row: [CA] = 0.1%

figure(556)
ind = 1;
for i = 2:6
    
    xValues = transpose(repmat(log10(xaxisValues),6,1));
    yValues = ([squeeze(simp_1(:,i,:)) squeeze(simp_2(:,i,:))]);

    subplot(1,5,7-i)
    hold on

    plot(xaxisValues,yValues,'.','color',[1,1,1]*0.7,'markersize',10)
    totSeq = (squeeze(simp_1(:,i,:))+squeeze(simp_2(:,i,:)))/2;

    plot(xaxisValues,mean(totSeq,2),'-','linewidth',1,'color',0.4*[1,1,1])
    errorbar(xaxisValues,mean(totSeq,2),std(totSeq,[],2),'k.','linewidth',1,'color',0.4*[1,1,1])
    plot(xaxisValues,mean(totSeq,2),'ko','linewidth',2,'markersize',8)
   
% MODEL THE TRENDS TO FIND STATISTICAL SIGNIFICANCE
    if i == 4 % correspond to the biphasic response; splitting the fit process
        % into 3 segments 
        
        modelX{i-1} = cell(2,1);
        xValuesLeft = xValues(1:2,:);
        yValuesLeft = yValues(1:2,:);
        xValuesMiddle = xValues(2:4,:);
        yValuesMiddle = yValues(2:4,:); 
        xValuesRight = xValues(4:5,:);
        yValuesRight = yValues(4:5,:);   
        %Error 
        modelX{i-1}{1} = fitlm(xValuesLeft(:),yValuesLeft(:));
        modelX{i-1}{2} = fitlm(xValuesMiddle(:),yValuesMiddle(:));        
        modelX{i-1}{3} = fitlm(xValuesRight(:),yValuesRight(:));
        
        slopeX(3,2) = modelX{i-1}{1}.Coefficients.pValue(2)/2;
        slopeX(3,1) = modelX{i-1}{1}.Coefficients.Estimate(2);
        slopeX(4,2) = modelX{i-1}{2}.Coefficients.pValue(2)/2;
        slopeX(4,1) = modelX{i-1}{2}.Coefficients.Estimate(2);        
        slopeX(5,2) = modelX{i-1}{3}.Coefficients.pValue(2)/2;
        slopeX(5,1) = modelX{i-1}{3}.Coefficients.Estimate(2);
        
        Xnew = xValues;

        ind = ind+3;
             
    else
        
        Xnew = xValues;        
        modelX{i-1} = fitlm(xValues(:),yValues(:));
        slopeX(ind,2) = modelX{i-1}.Coefficients.pValue(2)/2;
        slopeX(ind,1) = modelX{i-1}.Coefficients.Estimate(2);

        ind = ind+1;
    end
    
    hold on
    set(gca,'xscale','log','xticklabel','','yticklabel','')

    axis([1 5000 3 20])
    title(titles(i-1))
end

subplot(1,5,1)
set(gca,'xscale','log','xtick',[1 10 100 1000],'xticklabel',{'10^0','10^1','10^2','10^3'},'ytick',0:5:25,'yticklabel',0:5:25)
xlabel('# of partitions')
ylabel('Simpson index')

set(gcf,'position',[0 0 850 250])
slopeX


