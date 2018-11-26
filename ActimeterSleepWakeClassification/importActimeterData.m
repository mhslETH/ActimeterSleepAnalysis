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
% Description: Importing raw MHSL Actimeter data 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Importing raw MHSL Actimeter data 

%% Configuration

clear
% Add the path of your file in the line below.
filename = 'D:\Study_Data\ATP2-14\170913_1311_ATP2-14_ankle_ActivityTracker_23_EB90F4070988_acc.csv';

delimiter = ',';
startRow = 3;
fs = 12.5; %Hz
range = 4; %g
%% Format for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
datetime = dataArray{:, 1};
x = dataArray{:, 2};
y = dataArray{:, 3};
z = dataArray{:, 4};
date = char(datetime);
time = date(:,12:23);

x = x * 0.0005 * range;
y = y * 0.0005 * range;
z = z * 0.0005 * range;

time_sec = time(:,7:12);
time_sec = str2num(time_sec);
date = date(:,1:10);
time = cellstr(time);
date = cellstr(date);
Ns = size(time_sec);
Ns = Ns(1);

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;