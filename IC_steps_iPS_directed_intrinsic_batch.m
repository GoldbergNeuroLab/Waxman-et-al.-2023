
%% Analysis of SELECT intrinsic properties of iPS-directed neurons from IC steps

% Julie Merchant, Goldberg lab
% last updated 6.12.23
% designed for analysis for Elisa Waxman's 2D cortical directed diff protocol manuscript 

% Inputs: IC steps protocol, .abf files in current clamp in FOLDERS
    % E.g. 5 pA steps starting from -30 pA
    
% User input per file:
    % Fig 1 - click time bound w/in which to identify V_rest or V_hold (exclude spontaneous PSPs and/or spikes)
    
% Outputs: 
    % Single Excel file:
        % 1st sheet is Intrinsic Properties, each row is a file
            % Passive membrane properties
                % Vm (mV) (spontaneous or holding)
                % Holding current (pA)
                % Input resistance (mOhms)
            % Properties of single AP (rheobase spike)
                % Cutoff (mV) used to determine an AP
                % Rheobase sweep number
                % Rheobase (pA)
                % AP threshold (mV)
            % Properties of repetitive firing
                % Max steady state firing freq (Hz)
        % 2nd sheet is Spikes info, concatenated column pairs
            % per pair, 1st column is current steps, 2nd column is # spikes per current step

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
Results = zeros(length(fileInfo),8);
SpikeTable = []; 

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
    time(indices,1)=(indices./(1e6./si));  % divide datapoints by sampling rate (Hz) to get time in s for each datapoint
    time_ms = time.*1000;  % get time in ms for each datapoint

    % Create 2D matrices from full data matrix
    matrix_V = zeros(data_points,num_sweeps);
    matrix_I = matrix_V;

    for i = 1:num_sweeps
        matrix_V(:,i) = data(:,1,i);   %creates matrix of VOLTAGE (OUTPUT) (each sweep = column)
        matrix_I(:,i) = data(:,2,i);   %creates matrix of CURRENT (INPUT) (each sweep = column)
    end

    %% Identify V_rest or holding voltage, holding current, and current steps

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
    V_hold = mean(mean(matrix_V(range_start:range_end,1:5)));  %find avg voltage w/in time window in 1st 5 sweeps, then avg all of those
    I_hold = mean(mean(matrix_I(range_start:range_end,:)));  %find avg current w/in time window in each sweep, then avg all of those

    % Loop to identify current injection for each sweep
    I_steps = zeros(1,num_sweeps); %preallocate
    for i = 1:num_sweeps
        I_steps(i) = round(mean(matrix_I(I_on:I_off,i))-I_hold);
    end

    %% Identify indices and number of spikes in each sweep

    % Set threshold for defining an AP (value in mV it has to cross)
    cutoff = -10;  %may need to adjust

    spike_table = zeros(num_sweeps,3); %preallocate to store [Sweep] [Current (pA)] [# spikes]
    spike_i = {}; %create empty cell array to store indices of spikes for each sweep (to calc ISIs)

    % Loop through each sweep to get indices of each spike and populate spike_table
    for sweep = 1:num_sweeps
        v = matrix_V(I_on:I_off,sweep); %make variable for that sweep's voltage during the current injection
        k = find(v > cutoff); %make variable to find indices where voltage > cutoff
        if isempty(k) %if no part in sweep crosses cutoff, no indices recorded
            spike_is = [];
        else
            spike_is = [k(1)]; %#ok<*NBRAK> %if at least 1 part in sweep crosses cutoff, record the first index
            for w = 2:length(k) %for new variable representing all remaining indices where voltage > cutoff
                if k(w) ~= k(w-1)+1 %only record an index  that is NOT continuous with the previous (aka identifying breaks so you don't double count an AP)
                    spike_is = [spike_is, k(w)]; %#ok<*AGROW>
                end
            end
        end
        spike_i{sweep} = spike_is; %#ok<*SAGROW> %w/in loop, add new cell with indices of each spike in that loop 

        spike_table(sweep,1) = sweep; %sweep #
        spike_table(sweep,2) = I_steps(sweep); %current injection
        spike_table(sweep,3) = size(spike_i{sweep},2); %record # of columns (aka # of spikes) per sweep
    end

    %% Rheobase and # spikes at rheobase

    spike_sweeps = spike_table(spike_table(:,3)>0); %get all sweep #s that have at least 1 spike
    rheo_sweep = spike_sweeps(1); %get the first sweep (from above) that has at least 1 spike

    rheobase = spike_table(rheo_sweep,2); %get current injection (pA) during rheo_sweep

    spikes_at_rheobase = spike_table(rheo_sweep,3); %get # spikes during rheo_sweep

    %% Max steady-state firing frequency

    [max_spikes,max_sweep] = max(transpose(spike_table(:,3)));  %get sweep w/ greatest # spikes from spike table

    max_SSFF = max_spikes./(time(I_off)-time(I_on));  %in Hz

    %% Input resistance
    % calculated from 1st trace (assumed to be negative current injection)

    V1 = mean(matrix_V((find(time==(time(I_on)-0.05))):I_on,1));  %find mean voltage in 50 ms before current step
    V2 = mean(matrix_V(((I_on+I_off)/2):I_off,1));  %find mean voltage in 2nd half of current step

    I1 = I_hold;  %get mean current before current step
    I2 = mean(matrix_I(((I_on+I_off)/2):I_off,1));

    Rm = (V2-V1)./(I2-I1).*1000;  %find Rm in mOhms

    %% AP threshold

    % Isolate voltage data from I_on to peak of 1st AP at rheobase
    AP1_V = matrix_V(I_on:(I_on+spike_i{rheo_sweep}(1)),rheo_sweep);

    AP1_V_d = gradient(AP1_V)./(si*1e-3);  %get 1st derivative of AP1_V (mV/ms)

    % Method 1: where 1st derivative = 10 mV/ms
    AP_thresh1is = find(AP1_V_d(100:end) >= 10);  %get indices where dV/dt >= 10 mV/ms
    AP_thresh1i = AP_thresh1is(1)+100;  %Get index where dV/dt is first >= 10 mV/ms (+100 indices?)
    AP_thresh1 = AP1_V(AP_thresh1i);

    
    %% Format conditions, file names, and all results in matrices w/in loop
    
    Conditions(fn,1) = AllNames(Conditioni);
    FileNames(fn,1) = AllNames(FileNamei);
    Results(fn,:) = [V_hold, I_hold, Rm, cutoff, rheo_sweep, rheobase, AP_thresh1, max_SSFF];
    SpikeTable{fn} = [spike_table(:,2:3)];
  
    clf(1);
    
end

%% Format results into tables w/ column and row labels

% Format intrinsic properties results w/ column/row labels
Results_IntrinsicProps = array2table([Conditions FileNames Results],'VariableNames',{'Folder','File name','V_hold (mV)',...
    'I_hold (pA)','Rm (mOhms)','cutoff (mV)','Rheobase sweep #','Rheobase (pA)','AP threshold (mV)','max SSFF (Hz)'});

% Pad each spike_table 
    for i = 1:size(SpikeTable,2)
        sizes(i) = [size(SpikeTable{i},1)];
        padSize = max(sizes);
    end

    paddedSpikeTable = [];

    for i = 1:size(SpikeTable,2)
        padded = padarray(SpikeTable{i}, [padSize-size(SpikeTable{i},1)], NaN,'post');
        paddedSpikeTable = [paddedSpikeTable padded];
    end

% Format I-F results w/ column/row labels

SpikeColumns = [];
    for i = 1:fn
       column = [Conditions(i) FileNames(i)];
       SpikeColumns = [SpikeColumns column];
    end

Results_Spikes = array2table([SpikeColumns; paddedSpikeTable]);

%% Save Excel

filename = 'Results_Batch_IntrinsicProps_and_Spikes.xlsx';
writetable(Results_IntrinsicProps,filename,'Sheet','IntrinsicProps','WriteVariableNames',true);
writetable(Results_Spikes,filename,'Sheet','Spikes','WriteVariableNames',false);

%% Save workspace
save("Workspace_Batch.mat");
