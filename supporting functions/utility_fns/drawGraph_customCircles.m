function [pLocations,colorsEach] = drawGraph_customCircles(A_sim, delta_sim, abundanceVector,abundThres,thrReg,pLocations,colorsEach)
%DRAWGRAPH_CUSTOMCIRCLES This function draws nodes only using scatter plot
%and arrange the nodes in multiple concentric circles to enable compact
%visualization of a large number of nodes 
%   

n = size(A_sim,1);

gam11 = sum(-A_sim.*(A_sim<0),2); % negative
gam22 = sum(A_sim.*(A_sim>0),2); % positive

gam11_norm = (gam11-min(gam11))/(max(gam11)-min(gam11));
gam22_norm = (gam22-min(gam22))/(max(gam22)-min(gam22));

delta = 1-delta_sim; % convert delta to growth rate
delta_norm = (delta-min(delta))/(max(delta)-min(delta));

% numOfCircles = 5;
radiusOfCircle = 200;
% rads = linspace(200/numOfCircles,radiusOfCircle,numOfCircles);
circleRes = 200; 

colorsType = [0,0,0;
    0,0,1;
    0,1,0;
    0,1,1;
    1,0,0;
    1,0,1;
    1,1,0;
    0.9,0.9,0.9];

% % k = rads*2*pi/(circleRes/numOfCircles); % (circleRes/numOfCircles) is the interval
% across two adjacent circles

% figure(152)
% for i = 1:numOfCircles
%     numOfDots = 8*i;
%     polarscatter(linspace(0,2*pi,numOfDots+1)+pi/8,rads(i).*ones(1,numOfDots+1))
%     hold on
% %     polarplot(linspace(th(i)-pi/4,th(i),circleRes),radiusOfCircle.*ones(1,circleRes),'color',colorsType(i,:))
% 
% %     polarplot(th(i).*ones(2,1),[0,radiusOfCircle],'color',[0.7,0.7,0.7])
% end

% figure(133)

th = pi/4:pi/4:2*pi;

for i = 1:8
    polarplot(th(i).*ones(2,1),[0,radiusOfCircle],'color',colorsType(i,:))
    hold on
    polarplot(linspace(th(i)-pi/4,th(i),circleRes),radiusOfCircle.*ones(1,circleRes),'color',colorsType(i,:))

%     polarplot(th(i).*ones(2,1),[0,radiusOfCircle],'color',[0.7,0.7,0.7])
end


% thrReg = 0.10;

if nargin == 5      
    pLocations = zeros(2,n);
    colorsEach = zeros(n,3);
    for i = 1:n
        % the color [gam11_norm(i), gam22_norm(i), delta_norm(i)]
        catNum = gam11_norm(i)*4 + gam22_norm(i)*2 + delta_norm(i);
        tempAng = rand()*pi/4*0.8+th(catNum+1)-pi/4+pi/4*0.1; % angle
        tempR = sqrt(rand())*radiusOfCircle*0.8+radiusOfCircle*0.2; % radius
        
        distanceMetric = min(sum(abs(pLocations-repmat([tempAng;tempR],1,n)).*repmat([1/(2*pi); 1/radiusOfCircle],1,n),1));
        while distanceMetric < thrReg
            tempAng = rand()*pi/4*0.8+th(catNum+1)-pi/4+pi/4*0.1; % regenerate angle
            tempR = sqrt(rand())*radiusOfCircle*0.8+radiusOfCircle*0.2; % regenerate radius
            
            distanceMetric = min(sum(abs(pLocations-repmat([tempAng;tempR],1,n)).*repmat([1/(2*pi); 1/radiusOfCircle],1,n),1));
        end
        pLocations(1,i) = tempAng;
        pLocations(2,i) = tempR;
        colorsEach(i,:) = [gam11_norm(i),gam22_norm(i),delta_norm(i)];
    end
end
        
indPres = abundanceVector>abundThres;
polarscatter(pLocations(1,indPres), pLocations(2,indPres),abundanceVector(indPres)*50,...
    colorsEach(indPres,:),'filled','MarkerEdgeColor',[0.7,0.7,0.7])

axis off





