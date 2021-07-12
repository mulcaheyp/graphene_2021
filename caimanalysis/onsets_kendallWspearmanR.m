%% opening data
clear all;clc
% open xlsx files containing the delF/Fo time series for individual
% pyramidal cells
rec1m1 = xlsread('');
rec1m1(:,1) = [];
rec2m1 = xlsread('');
rec2m1(:,1) = [];
mouse1cells = size(rec1m1,2);

CaimFs = ; % set the calcium imaging acquisition rate
%% normalizing data
% normalization factor for each cell
normfactor = rec1m1(1,:);
% normalize the ata
normrec1m1 = (rec1m1-normfactor)./normfactor;
normrec2m1 = (rec2m1-normfactor)./normfactor;

clear normfactor normrec1 normrec2m1
%% filterthedata
% Savitzky-Golay Filtering
normfiltrec1m1 = sgolayfilt(normrec1m1,2,7);
normfiltrec2m1 = sgolayfilt(normrec2m1,2,7);

clear normrec1m1 normrec2m1
%% ONSETS
% Define windows are seizure onsets. In this animal, there are three onsets
% captured in the two records
onset1 = normfiltrec1m1(45:65,:);
onset2 = normfiltrec2m1(160:180,:);
onset3 = normfiltrec2m1(2905:2925,:);

entirerecords = {normfiltrecm1 normfiltrec2m1 normfiltrec2m1};
segments = {onset1 onset2 onset3};

firstindices = [45 160 2905];
% Get the onsets
for i = 1:length(segments)
    getrecord = entirerecords{i};
    getonset = segments{i};
    firstindex = firstindices(i);
    % Onset thresholds for every cell
    thresholds = mean(getrecord(firstindex-25:firstindex,:)) + 5*std((getrecord(firstindex-25:firstindex,:)));
    z = 1;
    for n = 1:length(thresholds)
        a = min(find(getonset(:,n)>thresholds(n)));
        if ~isempty(a)
            storeonsets(z,:) = [a n];
            z = z + 1;
        end
    end
    % Store the onsets for the seizure
    allonsets{i} = storeonsets;
    clear storeonsets z a 
    % Calculate and store the total recruitment time
    totalrecruitmenttime(i) = (max(storeonsets)-min(storeonsets))/CaimFs;
end

%%
% Group all the onsets for a single animal
mouse = {allonsets{1} allonsets{2} allonsets{3}};
clear allonsets entirerecords segments

%% Find cells that have been recruited to multiple seizures
% rank the cells in every seizure
for i = 1:length(mouse)
    mouse{i}(:,1) = mouse{i}(:,1)-min(mouse{i}(:,1));
    z = unique(mouse{i}(:,1));
    for n = 1:length(z)
        mouse{i}(mouse{i}(:,1)==z(n),3) = n;
    end
end
% store whether across seizures, individual cells were recruited or not
cellscommon_m1 = [1:mouse1cells]';
for i = 1:length(cellscommon_m1)
    for z = 1:length(mouse)
        giveme = mouse{z}(:,2);
        g = sum(giveme==i);
        storage_mouse1(i,z) = g;
    end    
end
% get the cells recruited in all three seizures
usethesecellsm1 = (sum(storage_mouse1')'==length(mouse));
cellstouse_m1 = cellscommon_m1(usethesecellsm1);
% bring their recruitment ranks into a matrix for KendallWspearmanR
% calculation
for i = 1:length(mouse)
    forKendallW_m1(:,i) = mouse{i}(mouse{i}(:,2)==cellstouse_m1,3);
end

%% Kendall's W and Spearman's R
% calculate Kendall's W
W1 = KendallCoef(forKendallW_m1);

% Generate a null distribution of Kendall's W
w1_null = [];
for i = 1:10000
    w1_null = [randi(max(forKendallW_m1(:,1)),size(forKendallW_m1,1),1) randi(max(forKendallW_m1(:,2)),size(forKendallW_m1,1),1) randi(max(forKendallW_m1(:,3)),size(forKendallW_m1,1),1)];
    w1_null_store(i) = KendallCoef(w1_null);
end
% p value calculation
pW1 = sum(W1<w1_null_store)/10000;


% calculate spearman R values
[spearman1,p1] = corr(forKendallW_m1,'Type','Spearman');

