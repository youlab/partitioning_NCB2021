function BI = simpsonInd (yend)
    yend = yend(yend>1E-9); 
    % 1e-9 is the minimum bar to count the population as present
    
    if sum(yend) == 0
        BI = 0;
    else
        
    totDens = sum(yend);
%     totDens
    p = yend./totDens;
%     p
    BI = 1./(sum(p.*p));
    end
end

