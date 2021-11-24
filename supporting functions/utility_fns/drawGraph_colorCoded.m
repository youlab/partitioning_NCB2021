function [G,p]=drawGraph_colorCoded(A_sim, delta_sim, names, abundanceVector, abundanceThreshold, colorMap ,graphType,normalizeFlag)

% this function is directly adapted from drawGraph_noEdges_moreEfficient().
% It does not plot the interaction edges at all. 
% This function does not use a pre-defined colorMap but uses model
% parameters to generate specific colors for each node. 

% A = A_sim'; % transpose interaction matrix used for simulation to interaction matrix for graph generation
n = size(A_sim,1);

gam11 = sum(-A_sim.*(A_sim<0),2); % negative
gam22 = sum(A_sim.*(A_sim>0),2); % positive
delta = 1-delta_sim; % convert delta to growth rate

if normalizeFlag ==1
    gam11 = (gam11-min(gam11))/(max(gam11)-min(gam11));
    gam22 = (gam22-min(gam22))/(max(gam22)-min(gam22));
    delta = (delta-min(delta))/(max(delta)-min(delta));
else
%     gam11 = gam11/(max(gam11)-min(gam11));
%     gam11_norm = gam11;
end

s = [1];
t = [1];

EdgeTable = table([s' t'], ...
    'VariableNames',{'EndNodes'});

% index = {1,2,3,4,5}';
% Ncolor = {colorMapPaperLight(1,:),colorMapPaperLight(3,:),colorMapPaperLight(5,:),colorMapPaperLight(7,:),colorMapPaperLight(9,:)}';
% NodeTable = table(names,index,Ncolor,'VariableNames',{'Name','Index','Ncolor'});
NodeTable = table(names,'VariableNames',{'Name'});
G = digraph(EdgeTable,NodeTable,'omitselfloops');

% G = digraph(A,names,'omitselfloops');

% hold off
p = plot(G,'Layout',graphType);
% p = plot(G,'Layout',graphType,'LineWidth',1*abs(G.Edges.Weight),'ArrowSize',5*abs(G.Edges.Weight));
% layout(p,'force','UseGravity',true)
% layout(p,'circle')
p.EdgeColor = [1,1,1];
p.LineWidth = 0.001; 

p.ArrowPosition =0.7;

% do not show the names
labelnode(p,1:n,{''})

for i = 1:n
%     abundanceVector
%     i
%     ind = mod(i,size(colorMap,1))+1;
%     ind
    if abundanceVector(i)> abundanceThreshold
        highlight(p,i,'NodeColor',[gam11(i), gam22(i), delta(i)])
        highlight(p,i,'MarkerSize',abundanceVector(i)*30)
    else % placeholder for population that is not present
        highlight(p,i,'NodeColor',[1,1,1])
        highlight(p,i,'MarkerSize',1)
    end
%     else
%         highlight(p,i,'NodeColor',[1,1,1])
%         nb_pre = predecessors(G,i); % neighbors
%         nb_suc = successors(G,i);
%         highlight(p,i,nb_pre,'EdgeColor',[1,1,1])
%         highlight(p,i,nb_suc,'EdgeColor',[1,1,1])
%     end
%     
end

% indsN = (G.Edges.Weight<0);
% for i = 1:size(G.Edges,1)
%     if indsN(i)== 1
%         highlight(p,G.Edges.EndNodes(i,:),'EdgeColor','r')
%     end
% end

% for i = 1:n
%     if abundanceVector(i)< abundanceThreshold
% 
%         highlight(p,i,'NodeColor',[1,1,1])
% %         nb_pre = predecessors(G,i); % neighbors
% %         nb_suc = successors(G,i);
% % %         nb_pre
% %         highlight(p,nb_pre,i,'EdgeColor',[1,1,1])
% %         highlight(p,nb_pre,i,'LineStyle','--')
% %         highlight(p,nb_pre,i,'LineWidth',0.001,'ArrowSize',0.001)
% %         
% %         highlight(p,i,nb_suc,'EdgeColor',[1,1,1])
% %         highlight(p,i,nb_suc,'LineStyle','--')
% %         highlight(p,i,nb_suc,'LineWidth',0.001,'ArrowSize',0.001)
%     end
% end

axis square
axis off



