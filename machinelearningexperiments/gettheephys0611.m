function feats0611 = gettheephys0611(x,lowphase,highampphase,highampamp,F)
    % calculate the electrophysiological features of interest
    % mean and standard deviation of the high gamma amplitude
    feat1 = mean(highampamp);
    feat2 = std(highampamp);
    
    % LFP and HG power
    [pxx,freqs] = periodogram(x,[],[],F);
    lowfpower = trapz(pxx(freqs>4&freqs<30));
    highgammapower = trapz(pxx(freqs>80&freqs<150));
    totalpower = trapz(pxx);
    
    % PLV
    plv = abs(1/length(x)*sum(exp(1i*lowphase-highampphase)));
    
    feats0611 = [feat1 feat2 lowfpower/totalpower highgammapower/totalpower plv];
end