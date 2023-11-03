%% Detecting presence/absence of sEPSCs from iPS-derived neuron voltage clamp traces

% Julie Merchant, Goldberg lab
% last updated 6.28.23
% designed for analysis of y/n sEPSCs for Elisa Waxman's 2D cortical directed diff protocol manuscript

% Inputs: VC gap-free protocol, .abf files in FOLDERS
    % generally -70 mV hold, 2 min recording per cell
    
% Outputs: single Excel file, each column is a cell/recording
    % first row is condition (day of diff), 2nd row is file name
    % remaining rows (3-123) are 0/1 sEPSCs (>6 SD of noise) detected in each 1 s chunk of the recording
    
% Dependencies: 
    % abfload in search path folder
    % .abf file(s) to be plotted in search path

%% Import files

% Get the directory of .abf files
dir_n = uigetdir('Select the directory:');  
fileInfo = dir(fullfile(dir_n,'**','*.abf'));

% Preallocate matrices for results parameters
Conditions = strings();
FileNames = strings();
Results = zeros(120,length(fileInfo));

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
    matrix_I = matrix_V; %#ok<NASGU>

    matrix_V = data(:,2);  %creates matrix of VOLTAGE (INPUT) 
    matrix_I = data(:,1);  %creates matrix of CURRENT (OUTPUT)
    
    %% Detect sePSCs 

    chunks = 1:(length(matrix_I)/120):length(matrix_I);

    epscs = zeros(length(chunks),1);

    for i=1:length(chunks)-1

        flipped = -matrix_I(chunks(i):chunks(i+1));
        flipped_clean = rmoutliers(flipped);
        threshold = 6*std(flipped_clean);
        baseline = mean(flipped_clean);
        if max(flipped) > baseline + threshold
            epscs(i,1) = 1;
        end

    end

    %% Format conditions, file names, and results in single matrix w/in loop
    
    Conditions(1,fn) = AllNames(Conditioni);
    FileNames(1,fn) = AllNames(FileNamei);
    Results(:,fn) = epscs;
    
end

%% Format results into table w/ column labels

Results_EPSCs = [Conditions; FileNames; Results];
filename = 'Results_Batch_EPSCs.xlsx';
writematrix(Results_EPSCs,filename);

