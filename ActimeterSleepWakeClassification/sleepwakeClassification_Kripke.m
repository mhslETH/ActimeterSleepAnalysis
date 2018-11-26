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
% Version:     v1.1
% Project:     Sleep Tracking with Accelerometer-Based Wearable Actimeter
% Description: Kripke Sleep Scoring Algorithm 
%              Preprocessing: Bandpass Filter (1.0-4.0Hz) and thresholding.
%              Activity Counts derived by Digital Integration per epoch.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Kripke Sleep Scoring Algorithm

%% Configuration

xyz = abs(x) + abs(y) + abs(z); % combine the axes

epoch_length = 60;              % [s] standard 60
thr = 0.08;                     % threshold for acceleration
label = 30;                     % epoch interval for labeling the x-axes
hp = 1.0;                       % [Hz] cut-off frequency highpass filter
lp = 4.0;                       % [Hz] cut-off frequency lowpass filter

Scoring_Threshold = 0.009;      % 0.009 for wrist, 0.0037 for ankle
Wake_Threshold = 0.019;         % 0.019 for wrist, 0.0037 for ankle
Rescoring_Threshold = 0.8;      % 0.8
Nonwear_Threshold = 0.03;       % 0.02
join_sleep_periods = 10;        % if less than X minutes awake, sleep periods are combined.
err = 4;                        % allowed sampling rate error

%% Highpass Filter
xyz_hp(1) = xyz(1);
alpha_hp = 1 / (2 * pi * (1/fs) * hp + 1);

for i = 2:Ns
    xyz_hp(i) = alpha_hp * (xyz_hp(i-1) + xyz(i) - xyz(i-1));
end

%Lowpass Filter
xyz_bp(1) = xyz_hp(1);
alpha_lp = (2 * pi * (1/fs) * lp) / (2 * pi * (1/fs) * lp + 1);

for i = 2:Ns
    xyz_bp(i) = xyz_bp(i-1) + alpha_lp * (xyz_hp(i) - xyz_bp(i-1));
end

%% Function to extract the sampling frequency for each second

sampling_frequencies = zeros(1,floor(Ns/(fs-5)));
sec = floor(time_sec(1));
s = 1;
pos = 1;

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

%% Thresholding to eliminate remaining noise

xyz_thr = xyz_bp;

for i=1:Ns
    if xyz_bp(i) < thr && xyz_bp(i) > -thr
        xyz_thr(i) = 0;
    end
end

%% Digital Integration per Epoch
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


%% Kripke 2010 Scoring Algorithm

activity_count = digital_integration;
j = size(activity_count);
for i = 6:(j(2)-2)
    D(i) = 0.3 * (0.0138*activity_count(i-5) + 0.0224*activity_count(i-4) + 0.0236*activity_count(i-3) + 0.0316 * activity_count(i-2) + 0.0944 *activity_count(i-1) + 0.03*activity_count(i) + 0.0212*activity_count(i+1));
end
j = size(D);
for i = 1:j(2)
    if D(i) > Scoring_Threshold
        sleep(i) = 1;
    else
        sleep(i) = 0;
    end
end

%% Determine Bedtime and Wake-up times

status = 0;
sleep_count = 1;
for i = 1:epochs - 6 - 15*(60/epoch_length)
    if status == 0
        if  nnz(activity_count(i:i+15*(60/epoch_length))) <= 1 && D(i) <= Wake_Threshold && sleep(i) ~= 0.8
            sleep_start_pos(sleep_count) = i;
            status = 1;
        end
    else
        if  D(i) >= Wake_Threshold && (sum(D(i:i+12*(60/epoch_length)) > Wake_Threshold)) > 8
            wake_up_pos(sleep_count) = i;
            status = 0;
            sleep_count = sleep_count + 1;
        end
    end
end

% if the measured wake period at the end is too short to detect wake-up 
% set wake_up at the end of recording.
j = size(sleep_start_pos);
j = j(2);
i = size(wake_up_pos);
i = i(2);
if j > i
    wake_up_pos(1,sleep_count) = (epochs - 8);
    sleep_count = sleep_count + 1;
end

%% Postprocessing to detect phases of low activity which are not sleep

rescored_sleep_count = 1;
for i = 1:sleep_count-1
    if nnz(D(sleep_start_pos(i):wake_up_pos(i))) > Rescoring_Threshold * (wake_up_pos(i)-sleep_start_pos(i))
    else
        rescored_sleep_start_pos(rescored_sleep_count) = sleep_start_pos(i);
        rescored_wake_up_pos(rescored_sleep_count) = wake_up_pos(i);
        rescored_sleep_count = rescored_sleep_count + 1;
    end
end

%% Count Sleep time during sleep periods

sleep_period_min = zeros(rescored_sleep_count-1,1);
sleep_rescoring = ones(1,epochs);

for i = 1:rescored_sleep_count-1
    sleep_period_min(i,1) = (rescored_wake_up_pos(i) - rescored_sleep_start_pos(i)) - sum(sleep(rescored_sleep_start_pos(i):rescored_wake_up_pos(i)));
    sleep_rescoring(rescored_sleep_start_pos(i):(rescored_wake_up_pos(i)-1)) = sleep(rescored_sleep_start_pos(i):(rescored_wake_up_pos(i)-1));
end

sleep_period = sleep_period_min;

%% Criteria to detect non-wear of the device during a detected sleep period
% just detects if the device lays still. (not if its in a moving bag)

for i = 1:rescored_sleep_count-1
    if nnz(activity_count(rescored_sleep_start_pos(i)+2:rescored_wake_up_pos(i)-5)) < Nonwear_Threshold * (rescored_wake_up_pos(i)-rescored_sleep_start_pos(i))
        sleep_rescoring(rescored_sleep_start_pos(i):rescored_wake_up_pos(i)) = 0.8;    
    end
end

%% Reduce displayed times to main sleep period.
% if Sleep_end and sleep_start are appart less than X minutes, they are
% combined for diaplay in table.

% extract time labels from data
for i = 1:epochs
    time_label(i) = time(sum(sampling_frequencies(1:epoch_length*i)));
    date_label(i) = date(sum(sampling_frequencies(1:epoch_length*i)));
end
time_label = char(time_label);
time_label = time_label(:,1:5);
for i = 1:floor(epochs/label)
    xtime(i,:) = time_label((i*label),:);
end
% get time labels
xtime(2:(floor(epochs/label)+1),:) = xtime(1:(floor(epochs/label)),:);
xtime(1,:) = time_label(1,:);

num_sleep = size(sleep_period);
date_sleep_start = date_label(rescored_sleep_start_pos(1:num_sleep(1)));
date_sleep_start = transpose(date_sleep_start);
rescored_sleep_start_pos = rescored_sleep_start_pos(1:num_sleep(1));
disp_rescored_wake_up_pos = rescored_wake_up_pos;
disp_rescored_sleep_start_pos = rescored_sleep_start_pos;
disp_date_sleep_start = date_sleep_start;
k = num_sleep(1)-1:-1:1;

for j = 1:num_sleep(1)-1
    i = k(j);
    if abs(rescored_wake_up_pos(i) - rescored_sleep_start_pos(i+1)) <= join_sleep_periods * (60/epoch_length)
        if nnz(activity_count(rescored_sleep_start_pos(i)+2:rescored_wake_up_pos(i)-5)) > Nonwear_Threshold * (rescored_wake_up_pos(i)-rescored_sleep_start_pos(i)) 
            if nnz(activity_count(rescored_sleep_start_pos(i+1)+2:rescored_wake_up_pos(i+1)-5)) > Nonwear_Threshold * (rescored_wake_up_pos(i+1)-rescored_sleep_start_pos(i+1))
                disp_date_sleep_start(i) = [];   
                disp_rescored_wake_up_pos(i) = [];
                disp_rescored_sleep_start_pos(i+1) = [];
            end
        end
    end       
end

%% Arrange Table to display results

% get activity counts for visualization and active minutes during sleep
% period, and sleep period time
activity_count_visualization = zeros(1,epochs);
n = size(disp_rescored_sleep_start_pos);
n = n(2);
minutes_moving = zeros(n,1);
for i = 1:n
    wake_time(i,1) = sum(sleep(disp_rescored_sleep_start_pos(i):disp_rescored_wake_up_pos(i)-1));
    minutes_moving(i) = nnz(activity_count(disp_rescored_sleep_start_pos(i):disp_rescored_wake_up_pos(i)-1));
    sleep_period_time(i,1) = (disp_rescored_wake_up_pos(i) - disp_rescored_sleep_start_pos(i));
    for j = disp_rescored_sleep_start_pos(i):disp_rescored_wake_up_pos(i)-1
          activity_count_visualization(j) = activity_count(j);
    end
    total_sleep_time(i,1) = sleep_period_time(i) - wake_time(i);
    if sleep_rescoring(disp_rescored_sleep_start_pos(i)) == 0.8
        minutes_moving(i) = 0;
        sleep_period_time(i,1) = 0;
        wake_time(i,1) = 0;
        total_sleep_time(i,1) = 0;
    end        
end

% Arrange a Table with the results
sleep_start = time_label(disp_rescored_sleep_start_pos,:);
wake_up = time_label(disp_rescored_wake_up_pos,:);
if num_sleep(1) > 1
    sleep_overview = table(disp_date_sleep_start, sleep_start, wake_up, total_sleep_time, sleep_period_time, minutes_moving, wake_time, ...
    'VariableNames',{'Date', 'Sleep_Start', 'Sleep_End','Total_Sleep_Time', 'Assumed_Sleep_Time','Minutes_Moving','Wake_Time'});
    disp(sleep_overview)
else
    disp('Date Sleep Start')
    disp(date_sleep_start)
    disp('Sleep_Start')
    disp(sleep_start)
    disp('Wake-up')
    disp(wake_up)
    disp('Assumed Sleep Time')
    disp(sleep_period)
    disp('Sleep Period Time')
    disp(sleep_period_time)
    disp('Minutes Moving')
    disp(minutes_moving)
    disp('Wake Time')
    disp(wake_time)
end

%% Visualization

figure;
plot(sleep_rescoring);
set(gca,'xTick',1:label:epochs,'XTickLabel',xtime)
xtickangle(45);
set(gca,'yTick',0:0.2:1, 'YTickLabel', {'Sleep'; ' '; ' ';' ';'Device Off';'Wake'})
axis normal
xlabel('Time')
title('Sleep/Wake with Kripke Scoring Algorithm')
ylim([-0.1 1.1])
hold on
plot(activity_count_visualization/5);
hold off

clearvars sleep_period_min sleep_period sleep_start sleep_start_pos  wake_up wake_up_pos sleep_overwiev wake_time ...
    rescored_sleep_start_pos rescored_wake_up_pos sleep_rescoring num_sleep sleep_overview sleep_period_time wake_time minutes_moving ...
    disp_date_sleep_start _pos disp_sleep_time date_sleep_start rescored_sleep_count ...
    s sec status err i j tspan xyz_hp sleep_count time_label disp_minutes_moving k n label lp hp ...
    disp_rescored_sleep_start_pos disp_rescored_wake_up
    
