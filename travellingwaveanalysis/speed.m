clear all;clc
% import the .mat file containing the gradients of the planes representing
% individual seizure discharges
load('')

% for each seizure, calculate the speeds from the first and last 25% of seizure discharges
for i = 1:length(all_seizure_gradients)
    % get seizure
    get = all_seizure_gradients{i};
    % get the first and last 25% of gradients
    first25 = get(1:floor(0.25*length(get)),[2 3]);
    last25 = get(floor(0.75*length(get)):end,[2 3]);
    % calculate speeds
    first25 = 1./sqrt(first25(:,1).^2 + first25(:,2).^2);
    last25 = 1./sqrt(last25(:,1).^2 + last25(:,2).^2);
    % Reject speeds > 2 m/s
    first25(first25>2) = [];
    last25(last25>2) = [];
    % store the speeds from the first and last 25% of seizure discharges
    storage{i,1} = first25;
    storage{i,2} = last25;
    clear first25 last25
end