function [thespeeds thedirections] = givemetheplane(eeg,good_channels,datatype,F)
    % Select array configuration, see travelingwaveanalysis code for more
    % information
    if datatype==1
        all_channel_locs = {[0 0] [1500e-6 0] [1000e-6 0] [500e-6 0] [1500e-6 500e-6]...
        [1000e-6 500e-6] [500e-6 500e-6] [1500e-6 1000e-6] [1000e-6 1000e-6] [500e-6 1000e-6]...
        [1500e-6 1500e-6] [1000e-6 1500e-6] [500e-6 1500e-6] [0 1500e-6] [0 1000e-6] []};
    else
        all_channel_locs = {[0 500e-6] [0 0] [1500e-6 0] [1000e-6 0] [500e-6 0] [1500e-6 500e-6]...
    [1000e-6 500e-6] [500e-6 500e-6] [1500e-6 1000e-6] [1000e-6 1000e-6] [500e-6 1000e-6]...
    [1500e-6 1500e-6] [1000e-6 1500e-6] [500e-6 1500e-6] [0 1500e-6] [0 1000e-6]};
    end
    
    % filter in the LFP
    [b1,a1] = butter(2,[4 30]/(F/2),'bandpass');
    seiz_seg = filtfilt(b1,a1,eeg);
    % Set threshold
    lowamp = abs(seiz_seg);
    thresholds = mean(lowamp) + 2*std(lowamp);
    % Detect peaks on each channel 
    locations = {};
    no_spikes = [];
    for n = 1:length(good_channels)
        [pks,locs] = findpeaks(lowamp(:,n),'MinPeakHeight',thresholds(n),'MinPeakDistance',0.2*2000);
        locations{n} = locs;
        no_spikes(n) = length(locs);
        clear locs pks
    end

    clear toamp IMAGS REALS thresholds

    % find spikes that are common across channels
    spiketime  =[];
    % set limit on number of spikes that can be detected
    [nospikestoconsider,id] = max(no_spikes);
    % set limit on indices in which to detect spike
    spikeiterative = [locations{id}-300 locations{id}+300];
    % If a channel has a peak within the indices defined above, store the
    % index, if not, store a 'NaN'
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
        % find how many NaNs are in a putative spike
        spike_temp = seizurespikes(i,:);
        A = sum(isnan(spike_temp));
        if A >= 3
        else
            spikes_for_store(q,:) = spike_temp;
            q = q + 1;
        end
        clear A spike_temp
    end 

    % model hippocampal space and peak timing
    for i = 1:size(spikes_for_store,1)
        % Get times of the peak LFP amplitude
        get_spikes = spikes_for_store(i,:);
        get_spikes = get_spikes'/2000;
        % associate channel number and peak LFP amplitude time
        associations = [good_channels' get_spikes];
        associations(isnan(associations(:,2)),:) = [];
        % associate channel location and peak LFP amplitude time
        for n = 1:size(associations,1)
            for_plane(n,:) = [cell2mat(all_channel_locs(associations(n,1))) associations(n,2)];
        end
        storage_for_plane{i} = for_plane;
        clear get_spikes associations for_plane
    end

    % fit models, find direction and speed of the wave
    for i = 1:length(storage_for_plane)
        for_model = storage_for_plane{i};
        filtmd = fit(for_model(:,[1 2]),for_model(:,3),'poly11');
        xcoeff = filtmd.p10;
        ycoeff = filtmd.p01;
        gradients(i,:) = [mean(for_model(:,3)) xcoeff ycoeff];
        clear for_model filtmd xcoeff ycoeff
    end
    
    % store the gradient directions and speeds
    thedirections = [gradients(:,1) atan2(gradients(:,3),gradients(:,2))];
    thespeeds = [gradients(:,1) 1./(sqrt(gradients(:,2).^2 + gradients(:,3).^2))];
end