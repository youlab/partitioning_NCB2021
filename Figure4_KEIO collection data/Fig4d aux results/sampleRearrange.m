function BIs_ = sampleRearrange(BIs_pre)

    % 1st dimension: partitioning level
    % 2nd [ca]
    % 3rd repeats
    
    BIs = transpose(reshape(BIs_pre,12,8));

    BIs_ = zeros(5,6,3);
    BIs_(:,:,1) = BIs(1:5,1:6);

    BIs_(1:3,:,2) = BIs(1:3,7:12);
    BIs_(4,:,2) = BIs(5,7:12);
    BIs_(5,:,2) = BIs(7,1:6);

    BIs_(1,:,3) = BIs(6,1:6);
    BIs_(2,:,3) = BIs(6,7:12);
    BIs_(3,:,3) = BIs(4,7:12);
    BIs_(4,:,3) = BIs(8,1:6);
    BIs_(5,:,3) = BIs(7,7:12);
    
