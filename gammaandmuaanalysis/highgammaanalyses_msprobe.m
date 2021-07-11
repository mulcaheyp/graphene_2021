clear all;clc
%% Initialization
%% Load channel labels and annotations
% import the channel labels and annotations for a given multishank probe
% record
load('')
load('')
%% Import experimental data and downsample to 2 kHz
% select electrophysiological data from one experiment and import into the work
% space
[records,pathname] = uigetfile('.rhd','Multiselect','on');
exp_data = [];
norecords = size(records,2);
for i = 1:norecords
    name = char(records(i));
    modded_read_Intan_RHD2000_file(name,pathname);
    downsampdat = downsample(amplifier_data',10);
    if i == 1
        exp_data = downsampdat;
    else
        exp_data(end+1:end+length(downsampdat),:) = downsampdat;
    end
    clear amplifier_data downsampdat
end
clear name downsampdat i amplifier_channels aux_input_channels 
clear aux_input_data filename frequency_parameters norecords notes path pathname
clear records spike_triggers supply_voltage_channels supply_voltage_data t_amplifier
clear t_aux_input t_supply_voltage

% generate a time matrix for entire, downsampled experimental data
time_mat_exp = 1/2000*[1:length(exp_data)];



%% Filter to generate low (4-30 Hz) and high gamma (80-150 Hz) signals
% create LFP and HG filters
[low_num,low_denom] = butter(2,[4 30]/(2000/2),'bandpass');
[hg_num,hg_denom] = butter(2,[80 150]/(2000/2),'bandpass');

lowfreq = zeros(size(exp_data));
highfreq = zeros(size(exp_data));
% filter all channels
for i = 1:32
    lowfreq(:,i) = filtfilt(low_num,low_denom,exp_data(:,i));
    highfreq(:,i) = filtfilt(hg_num,hg_denom,exp_data(:,i));
end

clear low_num low_denom hg_num hg_denom
%% Get seizures
% use the annotations object to get the times of the seizures within the
% records
getseizuretimes = annotations.annLayer;
for i = 1:size(getseizuretimes,1)
    annot = string(char(getseizuretimes(i,2)));
    if annot == "seizure"
        seizuretimes(i,:) = (getseizuretimes{i,1});
    end
end
seizuretimes_in_indices = seizuretimes*2000;
clear seizuretimes getseizuretimes annot i 

%% Analysis of the Data

%% Extraction of Phase Locking Behavior of High Gamma Activity
% Go through each seizure, get the phase locking behavior of the peak gamma
% amplitude from each channel in each seizure
for i = 1:size(seizuretimes_in_indices,1)
    % extract the LFP and HG segments for seizure
    hg_segment = highfreq(seizuretimes_in_indices(i,1):seizuretimes_in_indices(i,2),:);
    low_segment = lowfreq(seizuretimes_in_indices(i,1):seizuretimes_in_indices(i,2),:);
    time_segment = 1/2000*[1:length(low_segment)]';
    % calculate the LFP phase and the HG amplitude
    tophase = hilbert(low_segment);
    IMAGS = imag(tophase);
    REALS = real(tophase);
    lowphase = atan2(IMAGS,REALS);
    clear tophase IMAGS REALS
    tophase = hilbert(hg_segment);
    IMAGS = imag(tophase);
    REALS = real(tophase);
    highamp = sqrt(IMAGS.^2 + REALS.^2);    
    clear tophase IMAGS REALS
    % set thresholds for each channel
    thresholds = mean(highamp) + 2.5*std(highamp);
       
    for n = 1:32
        % on each channel, detect the peak high gamma amplitudes
        [pks,locs] = findpeaks(highamp(:,n),'MinPeakHeight',thresholds(n));
        % store the time and the amplitude of the peak high gamma
        % amplitudes
        HG_Amps{i,n} = [time_segment(locs) pks];
        % store the LFP phase at a detected peak high gamma amplitude event
        HG_PHASELOCKING{i,n} = [time_segment(locs) lowphase(locs,n)];
        % extract a null (randomly selected) distribution of LFP phases of
        % equal size to the experimentally observed distribution
        PHASELOCKING_null{i,n} = lowphase(randi(length(lowphase),length(locs),1));
        clear pks locs
    end
    clear hg_segment low_segment lowphase thresholds
end

% compile into an analysis object
analysis_obj.phaselocking = HG_PHASELOCKING;
analysis_obj.nulls = PHASELOCKING_null;
analysis_obj.info = 'Rows correspond to seizures. Columns to channels. Column 1 of each matrix is the time relative to the time relative to the beginning of the seizure segment, and column 2 is the Low freq phase.';

%% Saving the analysis object
% further compilation and saving to a .mat file
HG_Analyses.analysis = analysis_obj;

savename = char(strcat(string(annotations.date),'_gammaanalysis.mat'));

save(savename,'HG_Analyses')