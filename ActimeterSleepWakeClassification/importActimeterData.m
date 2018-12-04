%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2018 Mobile Health Systems Lab, ETH Zurich
% All rights reserved.  The use and distribution is permitted under the 
% 3-Clause BSD License
% https://opensource.org/licenses/BSD-3-Clause
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
filename = '..\data\MHSLActimeter_WristRecording.csv';

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
