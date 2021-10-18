function subgroups = createPartitionedVector(eachCell,groupInd,numPartitions)

    p = randperm(length(eachCell));
    eachCell_permed = eachCell(p);
    numSpecies = max(eachCell); 
    
    subgroups = zeros(numPartitions,numSpecies); 
    
    for i = 1:numPartitions
    
        cellsInGroup = eachCell_permed(groupInd==i);
        for j = 1:numSpecies
           subgroups(i,j) = sum(cellsInGroup==j);
        end
    end