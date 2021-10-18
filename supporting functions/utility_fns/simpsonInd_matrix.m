function BIs = simpsonInd_matrix (yend)
    n = size(yend,1);
    BIs = zeros(n,1);
    for i = 1:n
        BIs(i) = simpsonInd(yend(i,:));
    end
end