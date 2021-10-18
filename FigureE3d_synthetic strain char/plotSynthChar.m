% Created by Feilun Wu, feilunwu@gmail.com
% last edit: 1/15/2021

% This script plot the data from synthetic communities
% Extended Figure 3

%% initialize environment

close all
clear all

addpath(genpath('../supporting functions'))
setFigDef

%% read in data
[data,text] = xlsread('synth_comm_data.xlsx','Data');

endInd = 143;
od = data(4:33,1:endInd);
time = data(2,1:endInd);

%%
% only 3 out of the 5 strains were used and plotted

od_reshaped = reshape(od,6,5,143);

figure(323456)
ind = 1;
for i = [4 3 1]
    for j = [1 4 5]
        
        subplot(3,3,ind)
        plot(time/3600,squeeze(od_reshaped(j,i,:)),'color',[1,1,1]*0.4)
        axis([0 24 0 0.6])
        if ind ~= 7
            set(gca,'xticklabel','','yticklabel','')
        end
        ind = ind+1;
        
    end
end

subplot(3,3,7)
xlabel('Time (hour)')
ylabel('OD')

set(gcf,'position',[0 0 300 300])


