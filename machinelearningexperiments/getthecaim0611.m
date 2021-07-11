function featsfeats = getthecaim0611(X)
    % calculate the calcium imaging features of interest
    % mean cellular mean delF/Fo
    featureone = mean(mean(X));
    
    % mean cellular delF/Fo LL
    ddx = diff(X);
    featuretwo5 = mean(sum(abs(ddx))));
    % mean cellular d/dt delF/Fo standard deviation
    featurefour = mean(std(ddx));
    
    % mean and std cell-cell correlation
    corrs = corr(X);
    corrs(corrs==1) = [];
    feat5 = mean(corrs);
    feat6 = std(corrs);
    
    featsfeats = [featureone featuretwo5 featurefour feat5 feat6];
end