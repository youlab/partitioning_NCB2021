% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 1/15/2021

% This script plots pairwise interaction experimental data
% Extended Figure3f,g
% Figure 3d,e

%%
close all
clear all

addpath(genpath('../supporting functions'))
setFigDef

paperColors = paperColor; 

%%
data12 = xlsread('CFU_data.xlsx','negative pairwise');

cfu1 = data12(1:4,2:5);
cfu2 = data12(9:12,2:5);

% round of samples (total number of BI to be calculated)
numSample = 16;

BIHouse12 = zeros(numSample,4);
cfuSampled12 = zeros(numSample,2,4);

for j = 1:4 % 4 partitioning levels
    for i = 1:4 % 4 cfu1 plates -- total
        for k = 1:4  %cfu2 plates 
            thisTot = cfu1(j,i);
            this2 = cfu2(j,k);
            this1 = thisTot - this2 ;
            this1 = this1*(this1>0);
            BIHouse12(4 * (i-1) + k,j) = 1/((this1/thisTot)^2+(this2/thisTot)^2);

            cfuSampled12(4 * (i-1) + k,1,j) = this1/thisTot;
            cfuSampled12(4 * (i-1) + k,2,j) = this2/thisTot;
        end

    end
end 

%% ===============plot===============

N_ = [6,24,96,384];
figure(1)
subplot(1,2,1)
hold on
errorbar(N_,mean(BIHouse12),std(BIHouse12),'.','linewidth',1,'color',0.4*[1,1,1])
plot(N_,mean(BIHouse12),'-','linewidth',1,'color',0.4*[1,1,1])

set(gca,'xscale','log','xtick',[10,100,1000])
axis([2,1000,1,2.1])


plot(N_,mean(BIHouse12),'ok','linewidth',2,'markersize',10)
xlabel('# of partitions')
ylabel('Simpson index')
axis square

figure(2)
subplot(1,2,1)
hold on
comp12 = transpose(squeeze(mean(cfuSampled12,1)))*100;
std12 = squeeze(std(cfuSampled12,1))*100;
ar12 = area(comp12,'facecolor','r');
ar12(1).FaceColor = paperColors(7,:);
ar12(2).FaceColor = paperColors(3,:);
xticklabels({'6','24','96','384'})
errorbar(1:4,comp12(:,1),std12(1,:),'k.','linewidth',1);
ytickformat( 'percentage');
axis([1 4 0 100])
axis square


%% ===============plot===============
N_ = [6,24,96,384];
figure(1)
subplot(1,2,2)
hold on
errorbar(N_,mean(BIHouse13),std(BIHouse13),'.','linewidth',1,'color',0.4*[1,1,1])
plot(N_,mean(BIHouse13),'-','linewidth',1,'color',0.4*[1,1,1])

set(gca,'xscale','log','xtick',[10,100,1000])
axis([2,1000,1,2.1])

plot(N_,mean(BIHouse13),'ok','markersize',10,'linewidth',2)
xticklabels({''})
yticklabels({''})
axis square

figure(2)
subplot(1,2,2)
hold on
comp13 = transpose(squeeze(mean(cfuSampled13,1)))*100;
std13 = squeeze(std(cfuSampled13,1))*100;
ar13 = area(comp13,'facecolor','r');
ar13(1).FaceColor = paperColors(1,:);
ar13(2).FaceColor = paperColors(3,:);

errorbar(1:4,comp13(:,1),std13(1,:),'k.','linewidth',1);
xticklabels({''})
yticklabels({''})

axis([1 4 0 100])
axis square

figure(2)
set(gcf,'position',[0 0 600 400])
figure(1)
set(gcf,'position',[0 0 600 400])

