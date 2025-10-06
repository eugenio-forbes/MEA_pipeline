%%% MEA Pipeline Step 2: Extracting Waveforms 
%%% from Maxwell Spike Times for Combinato Clustering
%%%
%%% For every individual channel in a recording, this function will
%%% extract spike waveforms based on the spike times stored in
%%% Maxwell data. Only negative spikes are extracted.
%%%
%%% Spike waveforms are then upsampled, aligned based on
%%% peak amplitude timing, and savedin HDF5 format in the
%%% same manner that is done by Combinato prior to clustering
%%% the spikes into groups.

function n02_extract_spikes_combinato_mw(varargin)
if isempty(varargin)
    root_directory = '/path/to/MEA_pipeline/parent_directory';
    subject = 'SC000';
    folder = 'yyyy-mm-dd_network';
else
    root_directory = varargin{1};  %%% (character array) Parent directory containing /MEA_pipeline and /MEA_database.
    subject = varargin{2};         %%% (character array) Subject code in 'SC000' format.
    folder = varargin{3};          %%% (character array) Folder name with data in 'yyyy-mm-dd_scan-name#' format. Raw recording in /folder/raw/data.raw.h5
end

%%%Declare directories and files
data_directory = fullfile(root_directory, 'MEA_database', subject, folder);
info_directory = fullfile(data_directory, 'info');
raw_directory = fullfile(data_directory, 'raw');
kilosort_directory = fullfile(data_directory, 'kilosort'); %%% To get the raw, recentered, resorted int16 recordings
combinato_directory = fullfile(data_directory, 'combinato_mw');
if isfolder(combinato_directory)
    rmdir(combinato_directory, 's');%%% Delete any previously saved data
end
mkdir(combinato_directory);

matrices_directory = strrep(combinato_directory, 'combinato', 'matrices');
if ~isfolder(matrices_directory)
    mkdir(matrices_directory)
end
listFile = fullfile(info_directory, 'recording_list.mat');
infoFile = fullfile(info_directory, 'session_info.mat');

%%% Get mwFile name
mwFile = dir(fullfile(raw_directory, '*.raw.h5'));
mwFile = fullfile(mwFile.folder, mwFile.name);

%%% Get recording list and session info
load(listFile, 'recording_list');
load(infoFile, 'session_info');
n_recordings = session_info.n_recordings;
n_total_channels = session_info.n_total_channels;
start_framenos = recording_list.start_frameno;
stop_framenos = recording_list.stop_frameno;
sampling_rates = recording_list.sampling_rate;

%%% Get total length of all files if concatenated to create dummy threshold variable for combinato
length_samples = stop_framenos(end)-start_framenos(1)+1;
length_ms = ceil(length_samples*1000)/sampling_rates(1);
combinato_threshold = [0;length_ms;25];

%%% Get start times and time steps in ms
start_times = zeros(length(start_framenos), 1);
for idx = 2:length(start_framenos)
    start_times(idx) = start_framenos(idx)-start_framenos(1)*(1000/sampling_rates(idx-1));
end

%%% Initialize cell arrays to hold spike waveforms and times for every
%%% possibly used channe; a matrix to track how many recordings had spikes
%%% for each channel (first column positive spikes, second column negative spikes);
%%% and a matrix to track which recording index had data for each channel.
all_waveforms = cell(n_total_channels, 1);
all_times = cell(n_total_channels, 1);
all_counts = zeros(n_total_channels, 1);
all_indices = false(n_total_channels, n_recordings);

%Initialize parallel pool if it hasn't been created
poolobj = gcp('nocreate');
if ~isempty(poolobj)
    delete(poolobj);
end
parpool(36)

for idx = 1:n_recordings    
    %%% Get recording info
    recording_name = recording_list(idx, :).recording_name{:};
    n_channels = recording_list(idx, :).n_channels;
    n_samples = recording_list(idx, :).n_samples;
    sampling_rate = sampling_rates(idx);
    start_time = start_times(idx);
    start_frameno = start_framenos(idx);
    filter_frequencies = recording_list(idx, :).filter_frequencies{:};
   
    recFolder = fullfile(kilosort_directory, recording_name);
    recFile = fullfile(recFolder, 'raw_recording.bin');
    sortingFile = fullfile(recFolder, 'sorted_mapping.mat');
    load(sortingFile, 'rowIDs', 'channels');
    
    %%% Read raw recordings
    fid = fopen(recFile, 'r');
    raw_recordings = fread(fid, [n_channels, n_samples]); fclose(fid);
    
    %%% Load mw spikes
    spikes = h5read(mwFile, ['/data_store/', recording_name, '/spikes']);
    framenos = [spikes.frameno] - start_frameno + 1;
    spike_channels = [spikes.channel]; clear spikes;
    
    %%% For some reason there are spike framenos that are lower than the
    %%% start frameno of its respective recording
    within_recording = framenos > 0 & framenos <= n_samples;
    framenos = framenos(within_recording);
    spike_channels = spike_channels(within_recording);
   
    %%% For each recorded channel, get respective spike times
    spike_times = arrayfun(@(x) (framenos(spike_channels == x)*(1000/sampling_rate)) + start_time , channels, 'UniformOutput', false);
    clear framenos spike_channels

    %Butterworth filter for extraction (default from MW)
    min_pass = 300; %Hz
    max_pass = 7000; %Hz
    order = 4;
    [butt_b, butt_a] = butter(order, [min_pass, max_pass]/(sampling_rate/2), 'bandpass');

    %%% Set parameters for spike extraction
    params = struct;
    params.indices_per_spike = round(2*(sampling_rate/1000)); % Equivalent to 2ms
    params.peak_index = round(params.indices_per_spike*0.3); % Peak of waveform at one third of the 2ms interval
    params.max_spike_duration = 1.5; % Spikes wider than 1.5 ms are thrown out in combinato
    params.pre_indices = params.peak_index-1;
    params.post_indices = params.indices_per_spike-params.peak_index;
    params.butt_b = butt_b;
    params.butt_a = butt_a;
    params.sampling_rate = sampling_rate;
    
    spike_waveforms = cell(n_channels, 1);
    
    parfor jdx = 1:n_channels
        channel_recording = raw_recordings(jdx, :);
        temp_times = spike_times{jdx};
        if ~isempty(temp_times)
            [temp_waveforms, temp_times] = get_channel_spikes(channel_recording, params, temp_times);
            if ~isempty(temp_waveforms)
                spike_waveforms{jdx} = single(temp_waveforms);
                spike_times{jdx} = temp_times;
            end
        end
    end
    clear raw_recordings
    
    for jdx = 1:n_channels
        if ~isempty(spike_waveforms{jdx})
            rowID = rowIDs(jdx);
            all_waveforms{rowID} = [all_waveforms{rowID}, spike_waveforms{jdx}];
            all_times{rowID} = [all_times{rowID};spike_times{jdx}];
            all_counts(rowID) = all_counts(rowID)+1;
        end
    end
    %%% To track recording indices for every channel
    all_indices(rowIDs, idx) = true;
end

%%% Create separate files for sorting negative spikes
do_sort_neg = find(all_counts>0);
do_sort_neg = arrayfun(@(x) sprintf('MEA_%05d/data_MEA_%05d.h5', x, x), do_sort_neg, 'UniformOutput', false);
writecell(do_sort_neg, fullfile(combinato_directory, 'do_sort_neg.txt'));

%%% Only save to combinato formatted files used channels that did have
%%% detected spikes.
not_empty = all_counts > 0;
all_rowIDs = find(not_empty);
all_waveforms = all_waveforms(not_empty);
all_times = all_times(not_empty);

n_active_channels = length(all_rowIDs);

parfor idx = 1:n_active_channels
    rowID = all_rowIDs(idx);
    save_combinato_h5(all_waveforms{idx}, all_times{idx}, combinato_threshold, rowID, combinato_directory);
end
delete(gcp('nocreate'));

%%% Save information that tracks used channels with spikes, number of
%%% recordings with spikes and indices of recordings for each channel.
save(fullfile(matrices_directory, 'sorted_channels.mat'), 'all_rowIDs', 'all_counts', 'all_indices');

end

function [spike_waveforms, spike_times] = get_channel_spikes(channel_recording, params, spike_times)
spike_waveforms = [];
good_times = spike_times>(params.pre_indices+5) & spike_times<(length(channel_recording)-params.post_indices-5);
spike_times = spike_times(good_times);
if length(good_times) <= 3
    return %%% There were not enough good spike times.
end
clear channel_detect
channel_extract = filtfilt(params.butt_b, params.butt_a, double(channel_recording'))';
extract_indices = [spike_times-params.pre_indices-5, spike_times+params.post_indices+5];
spike_waveforms = zeros(params.indices_per_spike+10, size(extract_indices, 1));
for idx = 1:size(extract_indices, 1)
    spike_waveforms(:, idx) = channel_extract(extract_indices(idx, 1):extract_indices(idx, 2));
end
spike_waveforms = upsample_spikes(spike_waveforms, 3);
[spike_waveforms, ~] = align_spikes(spike_waveforms, (params.pre_indices+5)*3-2, 15, 15, 'neg');
% spikes = clean_spikes(spikes, params.peak_index);
spike_waveforms = downsample_spikes(spike_waveforms, 3, params.indices_per_spike);
end

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
%     max_index = find(spikes(:, idx)==min(spikes(:, idx)));
%     bad_indices(idx) = max_index ~= center;
% end
% cleaned_spikes = spikes(:, ~bad_indices);
% end

function downsampled_spikes = downsample_spikes(spikes, skip, num_points)
indices = 1:skip:num_points*skip;
downsampled_spikes = spikes(indices, :);
end

function save_combinato_h5(waveforms, times, threshold, rowID, combinato_directory)

artifacts = int8(zeros(length(times), 1));

combinatoFolder = fullfile(combinato_directory, sprintf('MEA_%05d', rowID));
if ~isfolder(combinatoFolder)
    mkdir(combinatoFolder);
end

combinatoFile = fullfile(combinatoFolder, sprintf('data_MEA_%05d.h5', rowID));
if isfile(combinatoFile)
    delete(combinatoFile);
end

h5create(combinatoFile, '/thr', size(threshold), 'Datatype', 'double');
h5write(combinatoFile, '/thr', threshold);

h5create(combinatoFile, '/pos/times', 1, 'Datatype', 'double');
h5create(combinatoFile, '/pos/spikes', [1 1], 'Datatype', 'single');
h5create(combinatoFile, '/pos/artifacts', 1, 'Datatype', 'int8');

if ~isempty(waveforms)
    h5create(combinatoFile, '/neg/times', length(times), 'Datatype', 'double');
    h5create(combinatoFile, '/neg/spikes', size(waveforms), 'Datatype', 'single');
    h5create(combinatoFile, '/neg/artifacts', length(artifacts), 'Datatype', 'int8');
    h5write(combinatoFile, '/neg/times', times);
    h5write(combinatoFile, '/neg/spikes', waveforms);
    h5write(combinatoFile, '/neg/artifacts', artifacts);
else
    h5create(combinatoFile, '/neg/times', 1, 'Datatype', 'double');
    h5create(combinatoFile, '/neg/spikes', [1 1], 'Datatype', 'single');
    h5create(combinatoFile, '/neg/artifacts', 1, 'Datatype', 'int8');
end

end