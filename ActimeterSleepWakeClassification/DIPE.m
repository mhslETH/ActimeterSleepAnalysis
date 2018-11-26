%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2018 Mobile Health Systems Lab, ETH Zurich
% All rights reserved.  The use and distribution is permitted under the 
% 3-Clause BSD License
% https://opensource.org/licenses/BSD-3-Clause
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author:      Rafael Küttel, Mobile Health Systems Lab, ETH Zurich.
% Date:        2017
% Version:     v1.0
% Project:     Sleep Tracking with Accelerometer-Based Wearable Actimeter
% Description: 1st level Algorithm
%              Digital Integraion per Epoch (DIPE) with preprocessing, 
%              bandpass filter (1.17-4.3Hz) and thresholding.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Digital Integraion per Epoch (DIPE)

%% Configuration

% Add up the axes
xyz = abs(x) + abs(y) + abs(z);

sampling_frequencies = zeros(1,floor(Ns/(fs-5)));
sec = floor(time_sec(1));
err = 2; %allowed sampling rate error (just choose high enough)
s = 1;
pos = 1;
epoch_length = 60;

thr = 0.08;
label = 50;
alpha_hp = 0.629680;
alpha_lp = 0.683686;

n = size(xyz);
n = n(1);

%% Highpass Filter
xyz_hp(1) = xyz(1);

for i = 2:n
    xyz_hp(i) = alpha_hp * (xyz_hp(i-1) + xyz(i) - xyz(i-1));
end

%Lowpass Filter
xyz_bp(1) = xyz_hp(1);

for i = 2:n
    xyz_bp(i) = xyz_bp(i-1) + alpha_lp * (xyz_hp(i) - xyz_bp(i-1));
end

%% Function to extract the sampling frequency for each second

while pos+(fs*err) <= Ns - 2
    for i = pos:(pos+(floor(fs*err)))
        if time_sec(i) < (sec + 1) && time_sec(i) >= sec
            sampling_frequencies(s) = sampling_frequencies(s) + 1;
        end
    end
    sec = sec + 1;
    if sec == 60
        sec = 0;
    end
    
    pos = pos + sampling_frequencies(s);
    s = s + 1;
end

epochs = floor(find(sampling_frequencies,1,'last') / epoch_length);

%% Thresholding

xyz_thr = xyz_bp;

for i=1:Ns
    if xyz_bp(i) < thr && xyz_bp(i) > -thr
        xyz_thr(i) = 0;
    end
end

%% Partial Integration
%The tspan is updated every epoch to account for changing sampling
%frequencies of the sensor.

digital_integration = zeros(1,epochs);
xyz_thr = abs(xyz_thr);
integration = zeros(size(x));

for i = 1:epochs
    tspan = 1/ sampling_frequencies((i-1)*epoch_length + 1);
    for j = (sum(sampling_frequencies(1:(i-1)*epoch_length)))+1 : sum(sampling_frequencies(1:i*epoch_length))
        integration(j) = xyz_thr(j) * tspan + ((xyz_thr(j+1) - xyz_thr(j)) / 2) * tspan;
        digital_integration(i) = integration(j) + digital_integration(i);
    end
end

%% Visualization

for i = 1:epochs
    time_label(i) = time(sum(sampling_frequencies(1:epoch_length*i)));
    date_label(i) = date(sum(sampling_frequencies(1:epoch_length*i)));
end
time_label = char(time_label);
time_label = time_label(:,1:8);
for i = 1:floor(epochs/label)
    xtime(i,:) = time_label((i*label),:);
end
xtime(2:(floor(epochs/label)+1),:) = xtime(1:(floor(epochs/label)),:);
xtime(1,:) = time_label(1,:);

j = size(xyz);
j = j(1);
figure;
plot(xyz);
set(gca,'xTick',1:50000:j,'XTickLabel',time(1:50000:j))
axis normal
xtickangle(45);

figure;
plot(digital_integration);
set(gca,'xTick',1:label:epochs,'XTickLabel',xtime)
axis normal
xtickangle(45);
title('Integration over each Epoch');

