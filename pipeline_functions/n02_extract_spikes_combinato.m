%%% MEA Pipeline Step 2: Extracting Spikes for Combinato Clustering
%%%
%%% For every individual channel in a recording, this function will
%%% filter the signal, extract spike waveforms, time them and save
%%% them in HDF5 format in the same manner that is done by Combinato
%%% prior to clustering the spikes into groups.

function n02_extract_spikes_combinato(varargin)
if isempty(varargin)
    root_directory = '/path/to/MEA_pipeline/parent_directory';
    subject = 'SC000';
    folder = 'yyyy-mm-dd_network';
else
    root_directory = varargin{1};  %%% (character array) Parent directory containing /MEA_pipeline and /MEA_database.
    subject = varargin{2};         %%% (character array) Subject code in 'SC000' format.
    folder = varargin{3};          %%% (character array) Folder name with data in 'yyyy-mm-dd_scan-name#' format. Raw recording in /folder/raw/data.raw.h5
end

%%% Declare directories and files
data_directory = fullfile(root_directory, 'MEA_database', subject, folder);
info_directory = fullfile(data_directory, 'info');
kilosort_directory = fullfile(data_directory, 'kilosort'); %%% To get the raw, recentered, resorted int16 recordings

combinato_directory = fullfile(data_directory, 'combinato');
if isfolder(combinato_directory)
    rmdir(combinato_directory, 's'); %%% Delete any previously saved data
end
mkdir(combinato_directory);

matrices_directory = strrep(combinato_directory, 'combinato', 'matrices');
if ~isfolder(matrices_directory)
    mkdir(matrices_directory)
end

list_file = fullfile(info_directory, 'recording_list.mat');
info_file = fullfile(info_directory, 'session_info.mat');

%%% Get recording list and session info
load(list_file, 'recording_list');
load(info_file, 'session_info');
n_recordings = session_info.n_recordings;
n_total_channels = session_info.n_total_channels;
start_frame_numbers = recording_list.start_frame_number;
stop_frame_numbers = recording_list.stop_frame_number;
stop_frame_numbers = stop_frame_numbers - start_frame_numbers(1);
start_frame_numbers = start_frame_numbers - start_frame_numbers(1);
sampling_rates = recording_list.sampling_rate;

%%% Get total length of all files if concatenated to create dummy threshold variable for combinato
length_samples = stop_frame_numbers(end) + 1;
length_ms = ceil(length_samples * 1000) / sampling_rates(1);
combinato_threshold = [0; length_ms; 25];

%%% Get start times in ms
start_times = zeros(length(start_frame_numbers), 1);
for idx = 2:length(start_frame_numbers)
    start_times(idx) = start_frame_numbers(idx) * (1000 / sampling_rates(idx - 1));
end

%%% Initialize cell arrays to hold spike waveforms and times for every
%%% possibly used channel; a matrix to track how many recordings had spikes
%%% for each channel (first column positive spikes, second column negative spikes);
%%% and a matrix to track which recording index had data for each channel.
all_waveforms = cell(n_total_channels, 2);
all_times = cell(n_total_channels, 2);
all_counts = zeros(n_total_channels, 2);
all_indices = false(n_total_channels, n_recordings);

%%% Restart initialization of parallel pool in case another one is currently open
parallel_pool_handle = gcp('nocreate');
if ~isempty(parallel_pool_handle)
    delete(parallel_pool_handle);
end
parpool(36)

for idx = 1:n_recordings    
    
    %%% Get recording info
    recording_name = recording_list(idx, :).recording_name{:};
    n_channels = recording_list(idx, :).n_channels;
    n_samples = recording_list(idx, :).n_samples;
    sampling_rate = sampling_rates(idx);
    start_time = start_times(idx);
   
    recording_folder = fullfile(kilosort_directory, recording_name);
    recording_file = fullfile(recording_folder, 'raw_recording.bin');
    sorting_file = fullfile(recording_folder, 'sorted_mapping.mat');
    load(sorting_file, 'row_IDs');
    
    %%% Read raw recordings
    file_id = fopen(recording_file, 'r');
    raw_recordings = fread(file_id, [n_channels, n_samples]);
    fclose(file_id);
    
    %%% Butterworth filter for detection
    min_pass = 300; %Hz
    max_pass = 1000; %Hz
    order = 6;
    [butterworth_b_detect, butterworth_a_detect] = butter(order, [min_pass, max_pass] / (sampling_rate / 2), 'bandpass');
   
    %%% Butterworth filter for extraction
    min_pass = 300; %Hz
    max_pass = 3000; %Hz
    order = 6;
    [butterworth_b_extract, butterworth_a_extract] = butter(order, [min_pass, max_pass] / (sampling_rate / 2), 'bandpass');
    
    %%% Set parameters for spike extraction
    extraction_parameters = struct;
    extraction_parameters.indices_per_spike = round(2 * (sampling_rate / 1000));              %%% Equivalent to 2ms
    extraction_parameters.peak_index = round(extraction_parameters.indices_per_spike * 0.3);  %%% Peak of waveform at one third of the 2ms interval
    extraction_parameters.max_spike_duration = 1.5;                                           %%% Spikes wider than 1.5 ms are thrown out in combinato
    extraction_parameters.pre_indices = extraction_parameters.peak_index - 1;
    extraction_parameters.post_indices = extraction_parameters.indices_per_spike - extraction_parameters.peak_index;
    extraction_parameters.butterworth_b_detect = butterworth_b_detect;
    extraction_parameters.butterworth_a_detect = butterworth_a_detect;
    extraction_parameters.butterworth_b_extract = butterworth_b_extract;
    extraction_parameters.butterworth_a_extract = butterworth_a_extract;
    extraction_parameters.sampling_rate = sampling_rate;
    
    temp_positive_waveforms = cell(n_channels, 1);
    temp_positive_times = cell(n_channels, 1);
    temp_negative_waveforms = cell(n_channels, 1);
    temp_negative_times = cell(n_channels, 1);
    
    parfor jdx = 1:n_channels
    
        channel_recording = raw_recordings(jdx, :);
        
        [positive_spike_waveforms, positive_spike_times] = get_channel_spikes(channel_recording, extraction_parameters, 'pos');
        [negative_spike_waveforms, negative_spike_times] = get_channel_spikes(channel_recording, extraction_parameters, 'neg');
        
        if ~isempty(positive_spike_waveforms)
            temp_positive_waveforms{jdx} = single(positive_spike_waveforms);
            temp_positive_times{jdx} = (positive_spike_times * (1000 / sampling_rate)) + start_time;
        end
        
        if ~isempty(negative_spike_waveforms)
            temp_negative_waveforms{jdx} = single(negative_spike_waveforms);
            temp_negative_times{jdx} = (negative_spike_times * (1000 / sampling_rate)) + start_time;
        end
        
    end
    
    clear raw_recordings
    
    for jdx = 1:n_channels
    
        if ~isempty(temp_positive_waveforms{jdx})
            row_ID = row_IDs(jdx);
            all_waveforms{row_ID, 1} = [all_waveforms{row_ID, 1}, temp_positive_waveforms{jdx}];
            all_times{row_ID, 1} = [all_times{row_ID, 1}; temp_positive_times{jdx}];
            all_counts(row_ID, 1) = all_counts(row_ID, 1) + 1;
        end
        
        if ~isempty(temp_negative_waveforms{jdx})
            row_ID = row_IDs(jdx);
            all_waveforms{row_ID, 2} = [all_waveforms{row_ID, 2}, temp_negative_waveforms{jdx}];
            all_times{row_ID, 2} = [all_times{row_ID, 2}; temp_negative_times{jdx}];
            all_counts(row_ID, 2) = all_counts(row_ID, 2) + 1;
        end
        
    end
    
    %%% To track recording indices for every channel
    all_indices(row_IDs, idx) = true;
    
end

%%% Create separate files for sorting jobs of positive and negative spikes
do_sort_pos = find(all_counts(:, 1) > 0);
do_sort_pos = arrayfun(@(x) sprintf('MEA_%05d/data_MEA_%05d.h5', x, x), do_sort_pos, 'UniformOutput', false);
writecell(do_sort_pos, fullfile(combinato_directory, 'do_sort_pos.txt'));

do_sort_neg = find(all_counts(:, 2) > 0);
do_sort_neg = arrayfun(@(x) sprintf('MEA_%05d/data_MEA_%05d.h5', x, x), do_sort_neg, 'UniformOutput', false);
writecell(do_sort_neg, fullfile(combinato_directory, 'do_sort_neg.txt'));

%%% Only save used channels that did have detected spikes to combinato formatted files.
not_empty = all_counts(:, 1) + all_counts(:, 2) > 0;
all_row_IDs = find(not_empty);
all_waveforms = all_waveforms(not_empty, :);
all_times = all_times(not_empty, :);

n_active_channels = length(all_row_IDs);

parfor idx = 1:n_active_channels

    row_ID = all_row_IDs(idx);
    save_combinato_h5(all_waveforms(idx, :), all_times(idx, :), combinato_threshold, row_ID, combinato_directory);
    
end

delete(gcp('nocreate'));

%%% Save information that tracks used channels with spikes, number of
%%% recordings with spikes and indices of recordings for each channel.
save(fullfile(matrices_directory, 'sorted_channels.mat'), 'all_row_IDs', 'all_counts', 'all_indices');

end


%%% Function sets a threshold amplitude for timing potential spikes,
%%% filters signal for detection and extracts spikes waveforms, excludes
%%% waveforms based on timing and duration, upsamples the waveforms,
%%% aligns waveforms peak amplitude timing, downsamples the waveforms. 
function [spike_waveforms, spike_times] = get_channel_spikes(channel_recording, extraction_parameters, sign)

spike_waveforms = [];
spike_times = [];

time_step = 1000 / extraction_parameters.sampling_rate;

detection_signal = filtfilt(extraction_parameters.butterworth_b_detect, extraction_parameters.butterworth_a_detect, double(channel_recording'))';
threshold = 5 * (median(abs(detection_signal)) / .6745); %%% 5 times the estimated standard deviation (median / .6475).

if strcmp(sign, 'pos')
    pass_threshold = detection_signal > threshold;
else
    pass_threshold = detection_signal < (-1 * threshold);
end

if sum(pass_threshold) == 0
    return
end

waveform_borders = find(diff(pass_threshold) ~= 0);
first_border = find(diff(pass_threshold) == 1, 1, 'first');

waveform_borders = waveform_borders(find(waveform_borders == first_border):end);

if mod(length(waveform_borders), 2) > 0
    waveform_borders = waveform_borders(1:end - 1);
end

waveform_borders = reshape(waveform_borders, 2, []);

length_okay = diff(waveform_borders, 1, 1) * time_step <= extraction_parameters.max_spike_duration;

if sum(length_okay) <= 3
    return; %%% There were not enough good spikes
end

waveform_borders = waveform_borders(:, length_okay);

spike_times = zeros(size(waveform_borders, 2), 1);

for idx = 1:size(waveform_borders, 2)
    interval = detection_signal(waveform_borders(1, idx):waveform_borders(2, idx));
    if strcmp(sign, 'pos')
        spike_times(idx) = find(interval == max(interval), 1, 'first') + waveform_borders(1, idx) - 1;
    else
        spike_times(idx) = find(interval == min(interval), 1, 'first') + waveform_borders(1, idx) - 1;
    end
end

spike_times = spike_times(2:end);

good_times = spike_times > (extraction_parameters.pre_indices + 5) & spike_times < (length(channel_recording) - extraction_parameters.post_indices - 5);
spike_times = spike_times(good_times);

clear detection_signal

extraction_signal = filtfilt(extraction_parameters.butterworth_b_extract, extraction_parameters.butterworth_a_extract, double(channel_recording'))';

extraction_indices = [spike_times - extraction_parameters.pre_indices - 5, spike_times + extraction_parameters.post_indices + 5];

spike_waveforms = zeros(extraction_parameters.indices_per_spike + 10, size(extraction_indices, 1));

for idx = 1:size(extraction_indices, 1)
    spike_waveforms(:, idx) = extraction_signal(extraction_indices(idx, 1):extraction_indices(idx, 2));
end

spike_waveforms = upsample_spikes(spike_waveforms, 3);
[spike_waveforms, ~] = align_spikes(spike_waveforms, (extraction_parameters.pre_indices + 5) * 3 - 2, 15, 15, sign);
% spikes = clean_spikes(spikes, extraction_parameters.peak_index);
spike_waveforms = downsample_spikes(spike_waveforms, 3, extraction_parameters.indices_per_spike);

end


%%% Function to upsample spike waveforms by a factor interpolating with a spline
function upsampled_spikes = upsample_spikes(spikes, factor)

n_features = size(spikes, 1);
n_spikes = size(spikes, 2);

upsampled_n_features = (n_features - 1) * factor + 1;
axis = 1:factor:upsampled_n_features;
upsampled_axis = 1:upsampled_n_features;

upsampled_spikes = zeros(upsampled_n_features, n_spikes);

for idx = 1:n_spikes
    upsampled_spikes(:, idx) = interp1(axis, spikes(:, idx), upsampled_axis, 'spline');
end

end

%%% Function to align the timing of spike waveform peak amplitudes
function [aligned_spikes, index_maximum] = align_spikes(spikes, center, low, high, sign)

n_features = size(spikes, 1);
n_spikes = size(spikes, 2);

aligned_spikes = zeros(n_features-low-high, n_spikes);
index_maxima = zeros(n_spikes, 1);

for idx = 1:n_spikes
    
    if strcmp(sign, 'pos')
        [~, index_maxima(idx)] = max(spikes(center-low:center+high, idx));
        index_maxima(idx) = index_maxima(idx)+center-low-1;
    else
        [~, index_maxima(idx)] = min(spikes(center-low:center+high, idx));
        index_maxima(idx) = index_maxima(idx)+center-low-1;
    end
    
    aligned_spikes(:, idx) = spikes(index_maxima(idx)-center+low+1:index_maxima(idx)-center+n_features-high, idx);
    
end

index_maximum = center - low;

end


%%% Below function not used by combinato anymore
% function [cleaned_spikes, bad_indices] = clean_spikes(spikes, center)
% n_spikes = size(spikes, 2);
% bad_indices = false(n_spikes, 1);
% for idx = 1:n_spikes
%     max_index = find(spikes(:, idx) == min(spikes(:, idx)));
%     bad_indices(idx) = max_index ~= center;
% end
% cleaned_spikes = spikes(:, ~bad_indices);
% end

%%% Function to downsample spike waveforms by a factor, and return n_datapoints	
function downsampled_spikes = downsample_spikes(spikes, factor, n_datapoints)

indices = 1:factor:n_datapoints * factor;
downsampled_spikes = spikes(indices, :);

end

%%% Function to save spike times and waveforms in the format used by Combinato
function save_combinato_h5(waveforms, times, threshold, row_ID, combinato_directory)

positive_waveforms = waveforms{1};
negative_waveforms = waveforms{2};

positive_times = times{1};
negative_times = times{2};

positive_artifacts = int8(zeros(length(positive_times), 1));
negative_artifacts = int8(zeros(length(negative_times), 1));

combinato_folder = fullfile(combinato_directory, sprintf('MEA_%05d', row_ID));
if ~isfolder(combinato_folder)
    mkdir(combinato_folder);
end

combinato_file = fullfile(combinato_folder, sprintf('data_MEA_%05d.h5', row_ID));
if isfile(combinato_file)
    delete(combinato_file);
end

h5create(combinato_file, '/thr', size(threshold), 'Datatype', 'double');
h5write(combinato_file, '/thr', threshold);

if ~isempty(positive_waveforms)
    h5create(combinato_file, '/pos/times', length(positive_times), 'Datatype', 'double');
    h5create(combinato_file, '/pos/spikes', size(positive_waveforms), 'Datatype', 'single');
    h5create(combinato_file, '/pos/artifacts', length(positive_artifacts), 'Datatype', 'int8');
    h5write(combinato_file, '/pos/times', positive_times);
    h5write(combinato_file, '/pos/spikes', positive_waveforms);
    h5write(combinato_file, '/pos/artifacts', positive_artifacts);
else
    h5create(combinato_file, '/pos/times', 1, 'Datatype', 'double');
    h5create(combinato_file, '/pos/spikes', [1 1], 'Datatype', 'single');
    h5create(combinato_file, '/pos/artifacts', 1, 'Datatype', 'int8');
end
if ~isempty(negative_waveforms)
    h5create(combinato_file, '/neg/times', length(negative_times), 'Datatype', 'double');
    h5create(combinato_file, '/neg/spikes', size(negative_waveforms), 'Datatype', 'single');
    h5create(combinato_file, '/neg/artifacts', length(negative_artifacts), 'Datatype', 'int8');
    h5write(combinato_file, '/neg/times', negative_times);
    h5write(combinato_file, '/neg/spikes', negative_waveforms);
    h5write(combinato_file, '/neg/artifacts', negative_artifacts);
else
    h5create(combinato_file, '/neg/times', 1, 'Datatype', 'double');
    h5create(combinato_file, '/neg/spikes', [1 1], 'Datatype', 'single');
    h5create(combinato_file, '/neg/artifacts', 1, 'Datatype', 'int8');
end

end