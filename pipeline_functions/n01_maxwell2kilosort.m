%%% MEA Pipeline, Step 1: Converting Maxwell Biosystems Data to Kilosort Format
%%%
%%% Information about the structure of Maxwell Biosystems microelectrode array
%%% data structure below.
%%% Channel numbers and electrode numbers in MEA recordings are unrelated to their
%%% positions in the microelectrode array.
%%% This function will first gather information about the microelectrode array
%%% from the recording file so that a new channel mapping can be created.
%%% New channel mapping will associate x and y positions to a vector row (row_ID) so that
%%% results can be input into a vector that can be reshaped into a matrix of the same
%%% shape as the microelectrode array, for easier data manipulation and plotting.
%%% Microelectrode array recording data is saved in uint16 format. This function will
%%% multiply the data by the least significant bit so that the data is in microvolts,
%%% and it will demean the signal so that it can be saved in the int16 format expected by kilosort.
%%% Also saves a channel map in the format expected by kilosort.
%%% Additionally will plot sample signal and power spectral densities for all channels in recording.

function n01_maxwell_to_kilosort(varargin)
if isempty(varargin)
    root_directory = '/path/to/MEA_pipeline/parent_directory';
    subject = 'SC000';
    folder = 'yyyy-mm-dd_network';
    do_plot = false;
else                               %%% May edit input arguments above to run in editor or put them in this order in code:
    root_directory = varargin{1};  %%% (character array) Parent directory containing /MEA_pipeline and /MEA_database.
    subject = varargin{2};         %%% (character array) Subject code in 'SC000' format.
    folder = varargin{3};          %%% (character array) Folder name with data in 'yyyy-mm-dd_scan-name#' format. Raw recording in /folder/raw/data.raw.h5
    do_plot = varargin{4};         %%% (true/false)      Whether to make a plot of sample of raw recordings.
end

%%% Declare directories
data_directory = fullfile(root_directory, 'MEA_database', subject, folder);
raw_directory = fullfile(data_directory, 'raw');
info_directory = fullfile(data_directory, 'info');
if ~isfolder(info_directory)
    mkdir(info_directory);
end
kilosort_directory = fullfile(data_directory, 'kilosort');
if ~isfolder(kilosort_directory)
    mkdir(kilosort_directory);
end

%%% MaxLab file name
maxwell_file = fullfile(raw_directory, 'data.raw.h5');

%%% Load all of the data into a struct
maxwell_struct = n00_maxwell_file_to_struct(maxwell_file);

%%% How the data is structured:
%%% maxwell_struct:
%%% .hdf_version                  (string specifying version of hierarchical data file)
%%% .mxw_version                  (string specifying version of Maxlab)
%%% .version                      (string specifying version of device?)
%%% .bits                         (??? empty)
%%% .data_store:                  (where all recording info is stored)
%%%     .data(0000 - n-1):        (each individual recording has an ID starting from 0)
%%%         .recording_id         (i-1, int32)
%%%         .spikes:              (since we are extracting them again with other parameters, discarding these)
%%%             .frameno          (n_spikes x 1, int64, sample # of spike peaks)
%%%             .channel          (n_spikes x 1, int32, channel # of each spike)
%%%             .amplitude        (n_spikes x 1, single, spike peak amplitudes)
%%%         .start_time           (int64 of date of recording start in posixtime format)
%%%         .stop_time            (int64 of date of recording stop in posixtime format)
%%%         .well_id              (int32 with number of well, in this case only 0)
%%%         .groups.routed:
%%%             .channels         (n_channels x 1, uint16, channel number)
%%%             .frame_nos        (n_samples, uint64)
%%%             .raw              (n_samples x n_channels, uint16, raw recordings, switching to path to read)
%%%             .triggered        (???, int32, empty)
%%%         .settings:
%%%             .gain             (512, uncertain if it is gain used in hardware given that it is not used in Maxlab code)
%%%             .hpf              (default highpass filter frequency 300Hz, discarded)
%%%             .lsb              (least significant bit: value in microvolts of each increment of uint16 of raw signal)
%%%             .mapping          (locations in MEA in micrometers from origin of an electrode given channel/electrode number)
%%%                 .channel      (n_channels x 1, channel number, keeping in case it can be used for network scans)
%%%                 .electrode    (n_channels x 1, electrode number, keeping in case it can be used for network scans)
%%%                 .x            (n_channels x 1, x position in micrometers, used in this code to get row_ID)
%%%                 .y            (n_channels x 1, y position in micrometers, used in this code to get row_ID)
%%%             .sampling         (sampling rate in Hz)
%%%             .spike_threshold  (amplitude for spike detection, default 0 ?, discarded)
%%% .recordings                   (in .h5 file contains links to .data(id) corresponding to each recording, discarded)
%%% .wellplate:
%%%     .id                       (string specifying id of device)
%%%     .version                  (string specifying version of device)
%%%     .well000:
%%%        .id                    (int, well name - 0, discarded)
%%%        .name                  (string well name '1', discarded)
%%%        .(experiment_variable) 'Potassium', 'Carbachol', 'Stim', 'control', 'DIV', 'Slice'
%%% .wells:                       (links to data(id) corresponding to recordings of a given well, discarded)
%%%     .well000:
%%%         .events               Timing of stimulation events     

%%% Create MEA_xy_matrix for mapping electrode recordings
MEA_xy_matrix = n00_get_xy_matrix(root_directory, subject, folder);
MEA_total_channels = size(MEA_xy_matrix, 1);
MEA_pitch = min(diff(unique(MEA_xy_matrix(:, 1))));

%%% Determine number of recordings after deleting spike only recordings
recording_name = fields(maxwell_struct.data_store);
spikes_only = false(length(recording_name), 1);
for idx = 1:length(recording_name)
    if ~isfield(maxwell_struct.data_store.(recording_name{idx}).groups, 'routed')
        maxwell_struct.data_store = rmfield(maxwell_struct.data_store, recording_name{idx});
        spikes_only(idx) = true;
    end
end
recording_name(spikes_only) = [];
n_recordings = length(recording_name);

%%% Gather info for each recording and place in table
recording_id = NaN(n_recordings, 1);
well_id = NaN(n_recordings, 1);
start_time = cell(n_recordings, 1);
stop_time = cell(n_recordings, 1);
sampling_rate = NaN(n_recordings, 1);
n_channels = NaN(n_recordings, 1);
n_samples = NaN(n_recordings, 1);
start_frame_number = NaN(n_recordings, 1);
stop_frame_number = NaN(n_recordings, 1);
triggered = NaN(n_recordings, 1);
gain = NaN(n_recordings, 1);
least_significant_bit = NaN(n_recordings, 1);

for idx = 1:n_recordings
    recording_id(idx) = maxwell_struct.data_store.(recording_name{idx}).recording_id;
    well_id(idx) = maxwell_struct.data_store.(recording_name{idx}).well_id;
    temp_start = maxwell_struct.data_store.(recording_name{idx}).start_time;
    start_time{idx} = datestr(datetime(temp_start/1000, 'convertfrom', 'posixtime', 'TimeZone', 'local', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    temp_stop = maxwell_struct.data_store.(recording_name{idx}).stop_time;
    stop_time{idx} = datestr(datetime(temp_stop/1000, 'convertfrom', 'posixtime', 'TimeZone', 'local', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    sampling_rate(idx) = maxwell_struct.data_store.(recording_name{idx}).settings.sampling;
    n_channels(idx) = length(maxwell_struct.data_store.(recording_name{idx}).groups.routed.channels);
    n_samples(idx) = length(maxwell_struct.data_store.(recording_name{idx}).groups.routed.frame_nos);
    start_frame_number(idx) = maxwell_struct.data_store.(recording_name{idx}).groups.routed.frame_nos(1);
    stop_frame_number(idx) = maxwell_struct.data_store.(recording_name{idx}).groups.routed.frame_nos(end);
    triggered(idx) = maxwell_struct.data_store.(recording_name{idx}).groups.routed.triggered;
    gain(idx) = maxwell_struct.data_store.(recording_name{idx}).settings.gain;
    least_significant_bit(idx) = maxwell_struct.data_store.(recording_name{idx}).settings.lsb;
end

%%% 3 variables below just to make table with these names
subject = repelem({subject}, n_recordings, 1);
folder = repelem({folder}, n_recordings, 1);

recording_list = table(subject, folder, recording_name, recording_id, ...
    well_id, start_time, stop_time, sampling_rate, n_channels, n_samples, ...
    start_frame_number, stop_frame_number, triggered, gain, least_significant_bit);
    
save(fullfile(info_directory, 'recording_list.mat'), 'recording_list');

clear well_id start_time stop_time start_frame_number stop_frame_number triggered gain recording_list

subject = subject{1}; 
folder = folder{1};

%%% Loop through recordings
all_row_IDs = NaN(sum(n_channels), 1);
channel_count = 0;

for idx = 1:n_recordings    

    %%% Gather raw recordings
    raw_path = maxwell_struct.data_store.(recording_name{idx}).groups.routed.raw;
    raw_recordings = h5read(maxwell_file, raw_path);
    
    %%% Create separate folder for each recording
    recording_directory = fullfile(kilosort_directory, recording_name{idx});
    if ~isfolder(recording_directory)
        mkdir(recording_directory);
    end
    
    %%% Gather recording mapping data
    channels = maxwell_struct.data_store.(recording_name{idx}).settings.mapping.channel;
    electrodes = maxwell_struct.data_store.(recording_name{idx}).settings.mapping.electrode;
    x = maxwell_struct.data_store.(recording_name{idx}).settings.mapping.x;
    y = maxwell_struct.data_store.(recording_name{idx}).settings.mapping.y;
    
    %%% Determine row numbers based on x and y positions for each channel
    row_IDs = arrayfun(@(i) find(MEA_xy_matrix(:, 1) == x(i) & MEA_xy_matrix(:, 2) == y(i)), (1:n_channels(idx))');
    [row_IDs, sorted_indices] = sortrows(row_IDs, 'ascend');
    
    %%% Resort data and save electrode IDs
    channels = channels(sorted_indices);
    electrodes = electrodes(sorted_indices);
    x = x(sorted_indices);
    y = y(sorted_indices);
    
    raw_recordings = raw_recordings(:, sorted_indices);
    
    save(fullfile(recording_directory, sprintf('sorted_mapping.mat')), 'row_IDs', 'channels', 'electrodes', 'x', 'y', 'sorted_indices');
    
    %%% Create kilosort channel mapping for this recording
    chanMap = 1:n_channels(idx);
    chanMap0ind = chanMap - 1;
    xcoords = x;
    ycoords = y;
    kcoords = ones(n_channels(idx), 1);
    connected = true(n_channels(idx), 1);
    fs = sampling_rate(idx);
    
    save(fullfile(recording_directory, 'chanMap.mat'), 'chanMap', 'chanMap0ind', 'xcoords', 'ycoords', 'kcoords', 'connected', 'fs')
    
    %%% Convert uint16 raw recordings to double and multiply by least_significant_bit
    %%% Also multiply by a million for the signal to be in uV
    raw_recordings = double(raw_recordings)*least_significant_bit(idx)*1000000;
    
    %%% Demean the signal to remove positive offset of uint16
    channel_means = mean(raw_recordings, 1);
    raw_recordings = raw_recordings - channel_means;
    
    %%% Round raw recording and convert to kilosort format int16 (n_channels x n_samples), so transpose matrix
    raw_recordings = int16(round(raw_recordings))';
    
    %%% Write recording matrix to int16 binary file.
    file_id = fopen(fullfile(recording_directory, 'raw_recording.bin'), 'wb');
    fwrite(file_id, raw_recordings, 'int16');
    fclose(file_id);
    
    if do_plot
    
        %%% Plot raw recording signal and power (500 ms sample)
        plot_duration = 500; %ms
        random_portion_idx = round(size(raw_recordings, 2)*(0.3 + (0.4*rand)));
        random_portion = random_portion_idx:random_portion_idx + (sampling_rate * (plot_duration / 1000)) - 1;
        raw_recordings = raw_recordings(:, random_portion);
        n00_plot_signal(raw_recordings, root_directory, subject, folder, recording_name{idx}, 'raw', sampling_rate(idx));
        clear raw_recordings
        
    end
    
    %%% Gather all row IDs to determine if this is an activity scan or
    %%% network recording.
    all_row_IDs(channel_count + 1:channel_count + n_channels(idx)) = row_IDs;
    channel_count = channel_count + n_channels(idx);
    
end

unique_IDs = unique(all_row_IDs);

%%% Determine type of scan
if length(unique_IDs) > 1024
    recording_type = {'activity'};
else
    recording_type = {'network'};
end

%%% Sumarize session info
subject = {subject};
folder = {folder};
hdf_version = maxwell_struct.hdf_version;
maxlab_version = maxwell_struct.mxw_version;
version = maxwell_struct.version;
wellplate_id = maxwell_struct.wellplate.id;
wellplate_version = maxwell_struct.wellplate.version;
recording_ids = {recording_id};
n_total_channels = MEA_total_channels;
n_channels_used = length(unique_IDs);

%%% Experiment variables
experiment_variables = fields(maxwell_struct.wellplate.well000);

if ismember('Potassium', experiment_variables)
    potassium = maxwell_struct.wellplate.well000.Potassium;
else
    potassium = NaN;
end

if ismember('Carbachol', experiment_variables)
    carbachol = maxwell_struct.wellplate.well000.Carbachol;
else
    carbachol = NaN;
end

if ismember('Stim', experiment_variables)
    if contains(maxwell_struct.wellplate.well000.Stim{:}, 'no')
        stim = true;
    else
        stim = false;
    end
else
    stim = false;
end

if ismember('control', experiment_variables)
    control = maxwell_struct.wellplate.well000.control;
else
    control = NaN;
end

if ismember('DIV', experiment_variables)
    DIV = maxwell_struct.wellplate.well000.DIV;
else
    DIV = NaN;
end

if ismember('Slice', experiment_variables)
    slice = maxwell_struct.wellplate.well000.Slice;
else
    slice = {''};
end

events = fields(maxwell_struct.wells.well000);

if ~isempty(events) && contains(events, 'events')
    stim_events = maxwell_struct.wells.well000.events;
    save(fullfile(info_directory, 'stim_events.mat'), 'stim_events');
end

complete = false;

%%% Make table and save
session_info = table(complete, subject, folder, hdf_version, maxlab_version, version, ...
    wellplate_id, wellplate_version, MEA_pitch, recording_type, n_total_channels, ...
    n_channels_used, n_recordings, recording_ids, potassium, carbachol, stim, ...
    control, DIV, slice);
    
save(fullfile(info_directory, 'session_info.mat'), 'session_info');

end