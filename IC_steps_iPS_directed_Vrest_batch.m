
%% Analysis of Vrest of iPS-directed neurons from IC steps (no hold)

% Julie Merchant, Goldberg lab
% last updated 6.12.23
% designed for analysis for Elisa Waxman's 2D cortical directed diff protocol manuscript 

% Inputs: IC steps protocol (no hold), .abf files in current clamp in FOLDERS
    % E.g. 5 pA steps starting from -30 pA
    
% User input per file:
    % Fig 1 - click time bound w/in which to identify V_rest or V_hold (exclude spontaneous PSPs and/or spikes)
    
% Outputs: 
    % Single Excel file with folder/condition/filename info and spontaneous Vrest (mV) 

% Dependencies:
    % abfload.m and extrema.m in MATLAB search path
    % .abf file(s) to be analyzed in MATLAB search path 
    
% Adaptations from:
    % SL_IC_current_steps_analysis_20230110.m
    % JPM_IC_steps_analysis_05_21_22_ala_09_19_22.m 
    % JPM_VC_HEK_activation_v9_ala.m
    % VC_HEK_activation_single.m

%% Open figure window
% Move figure window to other moniitor and make bigger, prior to proceeding

close all
clearvars 

set(0,'DefaultFigureWindowStyle','docked');
figure(1);

%% Import files

% Get the directory of .abf files
dir_n = uigetdir('Select the directory:');  
fileInfo = dir(fullfile(dir_n,'**','*.abf'));

% Preallocate matrices for results parameters
Conditions = strings();
FileNames = strings();
Results = zeros(length(fileInfo),1);

% Loop through each file and analyze
for fn = 1:length(fileInfo)
    pathFileName = fullfile(fileInfo(fn).folder,fileInfo(fn).name);
    AllNames = split(pathFileName,"\");
    FileNamei = size(AllNames,1);
    Conditioni = FileNamei-1;

    % Load the data
    [data,si,h] = abfload(pathFileName);
        % data = data frame configured as <data pts per sweep> by <number of channels> by <number of sweeps>.
        % si = sampling interval in microseconds (us)
        %       Note: sampling rate (Hz) = 1e6 / si (us)
        % h = information on file (selected header parameters)

    % Extract dimensions of full data matrix & set up loaded data
    [data_points, num_chan, num_sweeps] = size(data);
        % data_points = number of samples in each sweep (i.e. recorded data points) 
        % num_chan = number of channels
        % num_sweeps = number of sweeps in that file
    
    % Create vector w/ all indices in file and convert to time points
    indices = 1:data_points;
    time(indices,1)=(indices./(1e6./si));  %#ok<*SAGROW> % divide datapoints by sampling rate (Hz) to get time in s for each datapoint
    time_ms = time.*1000;  % get time in ms for each datapoint

    % Create 2D matrices from full data matrix
    matrix_V = zeros(data_points,num_sweeps);
    matrix_I = matrix_V;

    for i = 1:num_sweeps
        matrix_V(:,i) = data(:,1,i);   %creates matrix of VOLTAGE (OUTPUT) (each sweep = column)
        matrix_I(:,i) = data(:,2,i);   %creates matrix of CURRENT (INPUT) (each sweep = column)
    end

    %% Identify V_rest 
    % Identify when current step goes on/off
        %assumes step time is same for each sweep, and that last sweep is positive current injection
        
    I_diff = diff(matrix_I(:,end)); %extract dI/dt for last sweep (last step to maximize the difference from holding current)
    I_on = find(I_diff == max(I_diff),1); %get index when current goes on (max of dI/dt)
    I_off = find(I_diff == min(I_diff),1); %get index when current goes off (min of dI/dt)

    % Plot all traces (voltage) before I_on
    figure(1)
    plot(matrix_V(:,1:5));
        xlim([1 I_on+500]);
        xticks(1000:2000:I_on);
        ylim([-100 60]);
        xlabel('Time (ms)');
        ylabel('Voltage (mV)');

    % Define x axis range w/in which to find V_rest or V_hold
    [xi,yi] = getpts;
    range_start = round(xi(1));
    range_end = round(xi(2));

    % Identify holding voltage and current
    V_rest = mean(mean(matrix_V(range_start:range_end,1:5)));  %find avg voltage w/in time window in 1st 5 sweeps, then avg all of those

    %% Format conditions, file names, and results in matrices w/in loop
    
    Conditions(fn,1) = AllNames(Conditioni);
    FileNames(fn,1) = AllNames(FileNamei);
    Results(fn,:) = V_rest;
    
    clf(1);
    
end

%% Format results into tables w/ column and row labels

% Format intrinsic properties results w/ column/row labels
Results_Vrest = array2table([Conditions FileNames Results],'VariableNames',{'Folder','File name','V_rest (mV)'});

%% Save Excel

filename = 'Results_Batch_Vrest.xlsx';
writetable(Results_Vrest,filename,'Sheet','Vrest','WriteVariableNames',true);

%% Save workspace

save("Workspace_Batch_Vrest.mat");
