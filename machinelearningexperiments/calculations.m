clear all;clc
% load in a concatenated electrophysiological and caim seizure from a .mat
% file
load('')
% name used for saving the features
savename = 'formodel.mat';
% indicate the array configuration
recordno = 2;
%% plane fitting for traveling wave 

% extract the speeds for the traveling wave throughout the seizure
[speeds directs] = givemetheplane(multimodal_data_set.ephysdata,multimodal_data_set.channelguide,recordno,2000);
directs(speeds(:,2)>2,:) = [];
speeds(speeds(:,2)>2,:) = [];

%%
% generate a time matrix
time_mat = 1/2000*[1:length(multimodal_data_set.ephysdata)];

% calculate how many sliding windows will be included
NumWinsfcn = @(xLen, fs, winLen, winDisp) floor((xLen-winLen*fs)/(winDisp*fs)+1);
winL_insec = 4; % bin length in sec
winO_insec = 3; % bin shift in sec
no_wins = NumWinsfcn(size(multimodal_data_set.ephysdata,1),2000,winL_insec,winO_insec);

% determine the step sizes for the sliding windows in indices
stepsizeep = winO_insec*2000;
binsizeep = winL_insec*2000;
stepsizeca = floor(winO_insec*multimodal_data_set.CaImFs);
binsizeca = floor(winL_insec*multimodal_data_set.CaImFs);

% select the ephys channel closest to imaging plane
usethisephyschannel = find(multimodal_data_set.channelguide==multimodal_data_set.usechannel);
goodchannelephys = multimodal_data_set.ephysdata(:,usethisephyschannel);

% generate LFP and HG filters, clculate the LFP phase and HG amplitude
[b,a] = butter(2,[4 30]/(2000/2),'bandpass');
[b1,a1] = butter(2,[80 150]/(2000/2),'bandpass');
lowfreq = filtfilt(b,a,goodchannelephys);
highfreq = filtfilt(b1,a1,goodchannelephys);
tolowfreqphase = hilbert(lowfreq);
tohighfreqamp = hilbert(highfreq);
clear lowfreq highfreq b a b1 a1
lowfreqphase = atan2(imag(tolowfreqphase),real(tolowfreqphase));
highfreqamp = sqrt(imag(tohighfreqamp).^2 + real(tohighfreqamp).^2);
clear tolowfreqphase tohighfreqamp
tophaseofhighfreqamp = hilbert(highfreqamp);
phasehighfreqamp = atan2(imag(tophaseofhighfreqamp),real(tophaseofhighfreqamp));
clear tophaseofhighfreqamp 

% calculate features in the sliding windows
for i = 1:no_wins
    % get ephys and caim indices for the sliding windows
    idx1e = (i-1)*stepsizeep + 1;
    idx2e = (i-1)*stepsizeep + binsizeep;
    idx1c = (i-1)*stepsizeca + 1;
    idx2c = (i-1)*stepsizeca + binsizeca;
    % extract single channel ephys and multicellular delF/Fo data in the
    % window
    relevanttimes = time_mat(1,idx1e:idx2e);
    ephyssegment = multimodal_data_set.ephysdata(idx1e:idx2e,usethisephyschannel);
    caimsegment = multimodal_data_set.individualcells(idx1c:idx2c,:);
    
    % get the speeds within the window and find the median
    getthespeeds = speeds(speeds(:,1)>=relevanttimes(1)&speeds(:,1)<=relevanttimes(end),2);
    storethespeeds = median(getthespeeds);
    if isempty(storethespeeds)
        speedstorage(i,1) = NaN;
    else
        speedstorage(i,1) = storethespeeds;
    end
    clear getthespeeds storethespeeds
    
    % extract LFP phase and HG amplitude time series in the window
    lphse = lowfreqphase(idx1e:idx2e);
    highamp = highfreqamp(idx1e:idx2e);
    highampphase = phasehighfreqamp(idx1e:idx2e);
    
    % calculate and store ephys and caim features
    ephysfeaturesforstorage = gettheephys0611(ephyssegment,lphse,highampphase,highamp,2000);
    caimfeaturesforstorage = getthecaim0611(caimsegment);
    storeephys(i,:) = ephysfeaturesforstorage;
    storecaim(i,:) = caimfeaturesforstorage;
    clear ephysfeaturesforstorage caimfeaturesforstorage
end

%% Compile and save
% keep features only from bins that have a median speed
where = ~isnan(speedstorage);
goodephysfeatures = storeephys(~isnan(speedstorage),:);
goodcaimfeatures = storecaim(~isnan(speedstorage),:);
speedstorage(isnan(speedstorage))=[];

% compile and save to a .mat object
formodelbuilding.ephys = goodephysfeatures;
formodelbuilding.caim = goodcaimfeatures;
formodelbuilding.output = speedstorage;

save(savename,'formodelbuilding');