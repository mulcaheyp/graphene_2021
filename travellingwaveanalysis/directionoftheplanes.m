clear all;clc
% import the .mat file containing the gradients of the planes representing
% individual seizure discharges
load('')
recname = all_seizure_gradients;
clear all_seizure_gradients

% for each seizure in the record, calculate the angle of the gradient
for i = 1:length(recname)
    recnamedirect{i} = [recname{i}(:,1) atan2(recname{i}(:,3),recname{i}(:,2))];
end
clear i

% plot the distributions of gradient angles as polar histograms
figure()
for i = 1:length(recname)
    subplot(1,length(recname),i)
    polarhistogram(recnamedirect{i}(:,2))     
end

