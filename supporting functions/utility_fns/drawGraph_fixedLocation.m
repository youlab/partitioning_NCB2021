function [G,p]=drawGraph_fixedLocation(A_sim, names, abundanceVector, abundanceThreshold, colorMap ,graphType)

A = A_sim'; % transpose interaction matrix used for simulation to interaction matrix for graph generation
n = size(A,1);

N = length(names);

G = digraph(A,names,'omitselfloops');

% hold off
p = plot(G,'Layout',graphType,'LineWidth',1*abs(G.Edges.Weight),'ArrowSize',5*abs(G.Edges.Weight));
% layout(p,'force','UseGravity',true)
% layout(p,'circle')
p.EdgeColor = [78,117,163]/255;

p.ArrowPosition =0.7;
labelnode(p,1:N,{''})

for i = 1:n
%     abundanceVector
%     i
    if abundanceVector(i)> 0
        highlight(p,i,'NodeColor',colorMap(i,:))
        highlight(p,i,'MarkerSize',abundanceVector(i)*30)
    else % placeholder for population that is not present
        highlight(p,i,'NodeColor',colorMap(i,:))
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

indsN = (G.Edges.Weight<0);
for i = 1:size(G.Edges,1)
    if indsN(i)== 1
        highlight(p,G.Edges.EndNodes(i,:),'EdgeColor','r')
    end
end

for i = 1:n
    if abundanceVector(i)< abundanceThreshold

        highlight(p,i,'NodeColor',[1,1,1])
        nb_pre = predecessors(G,i); % neighbors
        nb_suc = successors(G,i);
%         nb_pre
        highlight(p,nb_pre,i,'EdgeColor',[1,1,1])
        highlight(p,nb_pre,i,'LineStyle','--')
        highlight(p,nb_pre,i,'LineWidth',0.001,'ArrowSize',0.001)
        
        highlight(p,i,nb_suc,'EdgeColor',[1,1,1])
        highlight(p,i,nb_suc,'LineStyle','--')
        highlight(p,i,nb_suc,'LineWidth',0.001,'ArrowSize',0.001)
    end
end

axis square
axis off

