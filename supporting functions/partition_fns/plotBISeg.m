function [BIseg, pres, densTot, BIsimpson]  = plotBISeg(D_,Y_end,xmax)

% xmax is the total carrying capacity for one population

M = length(D_);
BIseg = zeros(1,M);
pres  = zeros(1,M);
densTot = zeros(1,M);
BIsimpson = zeros(1,M);



for i = 1:M
    
    % take out the ones that are too low in abundance
    xmaxthis = xmax/D_(i);
    tee = Y_end{i}.*(Y_end{i}>1/xmaxthis);
    
    % sum all local environments
    tee_sum = sum(tee,1);
    
    BIseg(i) = shannonInd(tee_sum); 
    BIsimpson(i) = simpsonInd(tee_sum); 
    pres(i)  = sum(tee_sum >= 1/xmaxthis);
    densTot(i) = sum(tee_sum)./D_(i);
end


% subplot(2,1,1)% plot(D_,BIseg)
% set(gca,'xscale','log')

% subplot(2,1,2)
% plot(D_,pres)
% set(gca,'xscale','log')
