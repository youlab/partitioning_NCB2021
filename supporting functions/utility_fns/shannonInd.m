function BI = shannonInd (yend)
    yend = yend(yend>1E-9);
    if length(yend) ==1
        BI = 0;
    else
    totDens = sum(yend);
%     totDens
    p = yend./totDens;
%     p
    BI = -sum(p(p>0).*log(p(p>0)));
    end

