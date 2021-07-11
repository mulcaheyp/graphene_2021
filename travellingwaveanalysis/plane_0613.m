clear all;clc
%% Initialize
% load annotations for the record of interest
load('')

% open the data, downsample to 2 kHz, and concatenate
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

% create a time matrix to index the electrophysiological data
time_mat_exp = 1/2000*[1:length(exp_data)];

% get the times of the seizures in the record
getseizuretimes = annotations.annLayer;

for i = 1:size(getseizuretimes,1)
    annot = string(char(getseizuretimes(i,2)));
    if annot == "seizure"
        seizuretimes(i,:) = (getseizuretimes{i,1});
    end
end

% Convert to indices
seizuretimes_in_indices = seizuretimes*2000;
seizuretimes_in_indices(seizuretimes_in_indices(:,1)==0,:) = [];

% Extract channels with high quality electrophysiological data
good_channels = annotations.goodchans;
good_exp_data = exp_data(:,good_channels);
clear exp_data

% Select Electrode Configuration: This depends on the orientation of the
% implant relative to the custom ZIF connector and is determined during the
% experiment. One of the three Gr MEA records with high quality
% electrophysiological data had configuration 1, and the other two had
% configuration 2

% Configuration 1
% all_channel_locs = {[0 0] [1500e-6 0] [1000e-6 0] [500e-6 0] [1500e-6 500e-6]...
%     [1000e-6 500e-6] [500e-6 500e-6] [1500e-6 1000e-6] [1000e-6 1000e-6] [500e-6 1000e-6]...
%     [1500e-6 1500e-6] [1000e-6 1500e-6] [500e-6 1500e-6] [0 1500e-6] [0 1000e-6] []};
% Configuration 2
 all_channel_locs = {[0 500e-6] [0 0] [1500e-6 0] [1000e-6 0] [500e-6 0] [1500e-6 500e-6]...
    [1000e-6 500e-6] [500e-6 500e-6] [1500e-6 1000e-6] [1000e-6 1000e-6] [500e-6 1000e-6]...
    [1500e-6 1500e-6] [1000e-6 1500e-6] [500e-6 1500e-6] [0 1500e-6] [0 1000e-6]};
   
% Create LFP filter and filter the data
[b1,a1] = butter(2,[4 30]/(2000/2),'bandpass');
filted_good_exp_data = filtfilt(b1,a1,good_exp_data);


%% detect spikes
for zzzzzz = 1:length(seizuretimes_in_indices)
    % extract data from one seizure
    idxe1 = seizuretimes_in_indices(zzzzzz,1);
    idxe2 = seizuretimes_in_indices(zzzzzz,2);
    seiz_seg = filted_good_exp_data(idxe1:idxe2,:);
    
    % set threshold for detecting spikes
    lowamp = abs(seiz_seg);
    thresholds = mean(lowamp) + 2*std(lowamp);
    
    % detect spikes in the LFP on all channels with high quality
    % electrophysiological data
    locations = {};
    no_spikes = [];
    for n = 1:length(good_channels)
        [pks,locs] = findpeaks(lowamp(:,n),'MinPeakHeight',thresholds(n),'MinPeakDistance',0.2*2000);
        locations{n} = locs;
        no_spikes(n) = length(locs);
        clear locs pks
    end
    clear toamp IMAGS REALS thresholds

    % find spikes that occur across channels on the Gr MEA
    spiketime  =[];
    % set limit on number of spikes to expect to detect
    [nospikestoconsider,id] = max(no_spikes);
    % set limit on candidate indices for spikes occuring across the array
    spikeiterative = [locations{id}-300 locations{id}+300];
    % If a channel has a spike detected within the limits defined above,
    % store it, if a channel does not have a spike detected, store a 'NaN'
    % value
    for n = 1:length(spikeiterative)
        bounds = spikeiterative(n,:);
        for z = 1:length(good_channels)
            answer = locations{z}(find(locations{z}>=bounds(1,1) & locations{z}<=bounds(1,2)));
            if length(answer) > 1
                lmk = lowamp(answer,z);
                [~,idx] = max(lmk);
                answer = answer(idx);
            end
            if isempty(answer)
                spiketime(n,z) = NaN;
            elseif ~isempty(answer)
                spiketime(n,z) = answer;
            end
        end
    end
    seizurespikes = spiketime;

    % reject spikes that are marked with 'NaN' for too many channels
    q = 1;
    for i = 1:size(seizurespikes,1)
        % determine how many 'NaN's occur for a given spike
        spike_temp = seizurespikes(i,:);
        A = sum(isnan(spike_temp));
        if A >= 3
        else
            % If enough channels have a detected spikes, store for further
            % analysis
            spikes_for_store(q,:) = spike_temp;
            q = q + 1;
        end
        clear A spike_temp
    end 

    % model hippocampal space and peak timing
    for i = 1:size(spikes_for_store,1)
        % get the times of the peak amplitudes
        get_spikes = spikes_for_store(i,:);
        get_spikes = get_spikes'/2000;
        % associate the channels with the times of the peak amplitudes
        associations = [good_channels' get_spikes];
        associations(isnan(associations(:,2)),:) = [];
        % associate the channel locations with the times of the peak amplitudes
        for n = 1:size(associations,1)
            for_plane(n,:) = [cell2mat(all_channel_locs(associations(n,1))) associations(n,2)];
        end
        % store spatiotemporal information for fitting models
        storage_for_plane{i} = for_plane;
        clear get_spikes associations for_plane
    end

    % fit models, store the gradients of each of the planes
    for i = 1:length(storage_for_plane)
        for_model = storage_for_plane{i};
        filtmd = fit(for_model(:,[1 2]),for_model(:,3),'poly11');
        xcoeff = filtmd.p10;
        ycoeff = filtmd.p01;
        gradients(i,:) = [mean(for_model(:,3)) xcoeff ycoeff];
        clear for_model filtmd xcoeff ycoeff
    end
    
    % store the raw spike information, the plane models, and the gradients
    all_spikes_for_store{zzzzzz} = spikes_for_store;
    all_storage_for_plane{zzzzzz} = storage_for_plane;
    all_seizure_gradients{zzzzzz} = gradients;
    clear gradients clear storage_for_plane spikes_for_store
end

% Save the results (gradients) of the plane modeling to individual seizure discharges
save('','all_seizure_gradients')

