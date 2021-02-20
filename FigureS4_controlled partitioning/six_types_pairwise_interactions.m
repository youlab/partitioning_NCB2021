% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 1/15/2021

% This script simulates and plots the impact of partitioning on all 6 types
% of pairwise interaction
% We simplified the analysis by using two extremes: co-culture and
% mono-culture. 

%% initialize environment

close all
clear all

addpath(genpath('../supporting functions'))
setFigDef

%% initialize parameters

del_ = [0 1.2];
g_ =  [-2 0 10];

% % color scheme: 
paper_colors = paperColor();

colors =  [        
    0         0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

tend = 1000;
% set up the different segregation levels
D0 = 0.1;
xTot = 1e6;
N = 2;

options=odeset('NonNegative',1:N);
y0_orig = ones(1,N)/N;

rng(11)
[ y0_, typeCount ] = seedInit_comb(y0_orig,xTot,D0);

yend_house = cell(length(del_),length(del_),length(g_),length(g_));


for i = 1:length(del_)
    for j = 1:length(del_)
        for k = 1:length(g_)
            for l = 1:length(g_)
             
                del = [del_(i); del_(j)*0.9]; 
                % 0.9 is to avoid simulation artifacts when two parameters
                % are exactly the same
                gam = [0 g_(k); g_(l)*0.9 0];
                
                gam11 = gam.*(gam <0);
                gam22 = gam.*(gam >0);
                
                yend_house{i,j,k,l} = runSeg_2gamma({del, gam11, gam22}, y0_, typeCount, xTot, tend);

                % =================== plot the results ===================
                % identify which column the plot needs to go in
                if del_(i)>1 && del_(j)>1
                    col = 4;
                elseif (del_(i)>1 && del_(j)<1)
                    col = 3;
                elseif (del_(i)<1 && del_(j)>1)
                    col = 2;
                elseif (del_(i)<1 && del_(j)<1)
                    col = 1;
                else
                    col = NaN;
                end
                
                % identify which row the plot needs to go in
                if g_(k)==0 && g_(l)==0
                    rol = 1;
                elseif g_(l)>0 && g_(k)==0
                    rol = 2;
                elseif g_(l)<0 && g_(k)==0
                    rol = 3;
                elseif g_(l)>0 && g_(k)>0
                    rol = 4;
                elseif g_(l)>0 && g_(k)<0
                    rol = 5;
                elseif g_(l)<0 && g_(k)<0
                    rol = 6;
                else 
                    rol = NaN;
                end
                
               if ~isnan(rol)
                    [~,pres,~,BIsim] = plotBISeg(1:N,yend_house{i,j,k,l},xTot);

                    figure(8)
                    subplot(6,4,(rol-1)*4+col)
                    hold on

                    ptemp = bar([(yend_house{i,j,k,l}{1});sum(yend_house{i,j,k,l}{2})],'stacked');
                    set(ptemp(1),'FaceColor',paper_colors(9,:))
                    set(ptemp(2),'FaceColor',paper_colors(7,:))
                    set(gca,'yTick',[0 1 2],'xticklabel','','yticklabel','')
                    ptemp = plot([1,2],BIsim,'k-o');
                    
                    % depending on the trend, use different colors
                    if BIsim(1)> BIsim(2)
                        set(ptemp,'Color',colors(1,:));
                    elseif BIsim(1)== BIsim(2)
                        set(ptemp,'Color',[1,1,1]*0.7);
                    elseif BIsim(1)<= BIsim(2)
                        set(ptemp,'Color',colors(2,:));
                    end
                    axis([0.0 3.0 0 2.5])
               end
               % ========================================================

            end
        end
    end
end


subplot(6,4,21)
set(gca,'xtick',[1 2],'xticklabel',[1 2],'yticklabel',[0 1 2])
 