clear all;clc
%% Initialize
% load annotations and channel info in the record of interest
load('')
load('')

% select electrophysiological data from one experiment
[records,pathname] = uigetfile('.rhd','Multiselect','on');

% get seizure times in the record
getseizuretimes = annotations.annLayer;

for i = 1:size(getseizuretimes,1)
    annot = string(char(getseizuretimes(i,2)));
    if annot == "seizure"
        seizuretimes(i,:) = (getseizuretimes{i,1});
    end
end

% convert seizure times to indices
seizuretimes_in_indices = seizuretimes*20000;
clear seizuretimes getseizuretimes annot i 

% get the first index of every seizure
for n = 1:length(records)  
    beginningrecording_indices(n,1) = 1 + (n-1)*size(firsttenminutes,1);
end

clear n
%% Analysis and Write to External File
% create the spike band and LFP filters
[spikefiltnum,spikefiltdenom] = butter(2,[300 3000]/(20000/2),'bandpass');
[lownum,lowdenom] = butter(2,[4 30]/(20000/2),'bandpass');
% Through every seizure: import the data (not downsampled), filter, detect
% MUA, and extract the phase locking behavior on all functional MS probe channels
for n = 1:size(seizuretimes_in_indices,1)
    % determine which records to open to import the seizure
    seizureidx = seizuretimes_in_indices(n,:);
    A = seizureidx(1,1) >= beginningrecording_indices;
    record1open = max(find(A==1));
    clear A
    B = seizureidx(1,2) <= beginningrecording_indices;
    record2open = min(find(B==1));
    clear B
    if record1open==size(records,2)
        record2open=record1open;
    end
    if isempty(record2open)
        record2open = size(records,2);
    end
    % open the data and import into the MATLAB space
    i = 1;
    exp_data = [];
    for z = record1open:record2open
        name = char(records(z));
        modded_read_Intan_RHD2000_file(name,pathname);
        clear amplifier_channels aux_input_channels 
        clear aux_input_data filename frequency_parameters norecords notes path 
        clear spike_triggers supply_voltage_channels supply_voltage_data t_amplifier
        clear t_aux_input t_supply_voltage
        if i == 1
            exp_data = amplifier_data';
        else
            exp_data(end+1:end+length(amplifier_data),:) = amplifier_data';
        end
        clear amplifier_data 
        i = i + 1;
    end
    % extract the seizure from the opened data
    seizure_data = exp_data(seizureidx(1,1)-beginningrecording_indices(record1open):seizureidx(1,2)-beginningrecording_indices(record1open),:);
    clear exp_data record1open record2open
    % for all functional channels, detect MUA and extract phase locking
    % behavior
    for g = 1:32
        chan = seizure_data(:,g);
        % filter
        spikeband = filtfilt(spikefiltnum,spikefiltdenom,chan);
        % threshold and detect peaks
        muathresh = 4*median(abs(spikeband)/0.6745);
        spikeband = -spikeband;
        [pks,locs] = findpeaks(spikeband,'MinPeakHeight',muathresh);
        clear spikeband muathresh spikeband
        % calculate LFP phase
        low = filtfilt(lownum,lowdenom,chan);
        tophase = hilbert(low);
        IMAGS = imag(tophase);
        REALS = real(tophase);
        lowphase = atan2(IMAGS,REALS);
        clear tophase IMAGS REALS 
        % extract the index and LFP phase at detected MUA
        item2{n,g} = [locs lowphase(locs)];
        % extract a null phase locking distribution with equal size as the
        % experimentally observed MUA distribution
        nulldist{n,g} =  lowphase(randi(length(lowphase),length(locs),1));
        clear chan lowphase locs pks low
    end

end

%% Compile and save to a .mat file
analysis3_obj.muaphase = item2;
analysis3_obj.nulls = nulldist;
analysis3_obj.Note = 'Rows are seizures. Columns are channels. Col 1 of each cell is time of MUA, Col2 is Phase locking 4-30 Hz.';

MUA_Analyses.analysis3 = analysis3_obj;


savename = char(strcat(string(annotations.date),'_MUAanalysis.mat'));

save(savename,'MUA_Analyses')


