%%% MEA Pipeline, Step 6: Gather Kilosort Results
%%%
%%% In development. Function meant for spike extraction
%%% based on Phy processing of Kilosort results. 

function n06_get_kilosort_units(varargin)
if isempty(varargin)
    root_directory = '/path/to/MEA_pipeline/parent_directory';
    subject = 'SC000';
    folder = 'yyyy-mm-dd_network';
else
    root_directory = varargin{1};  %%% (character array) Parent directory containing /MEA_pipeline and /MEA_database.
    subject = varargin{2};         %%% (character array) Subject code in 'SC000' format.
    folder = varargin{3};          %%% (character array) Folder name with data in 'yyyy-mm-dd_scan-name#' format. Raw recording in /folder/raw/data.raw.h5
end

%%% Modifiable parameters

%%% Templates file is n_templates x n_samples x n_channels in size. For a
%%% given template ID, phy gets the templates pertaining to that ID for
%%% every channel. It then calculates the amplitude for each channel
%%% (maximum-minimum of the templates). The best channel for a template ID
%%% is the channel with the largest amplitude. In addition to this channel
%%% it selects, n_best_channels -1 closest channels to this one, if these
%%% other channels had signal, and the amplitude of their template exceeds
%%% a given threshold (size relative to the maximum template).
n_best_channels = 12; %%% Default: 12
amplitude_threshold = 0; %%% Default: 0

%%% Declare directories
data_directory = fullfile(root_directory, 'MEA_database', subject, folder);
info_directory = fullfile(data_directory, 'info');
kilosort_directory = fullfile(data_directory, 'kilosort'); %%% To get the raw, recentered, resorted int16 recordings
unit_directory = fullfile(data_directory, 'kilosort_units');
if ~isfolder(unit_directory)
    mkdir(unit_directory);
end

list_file = fullfile(info_directory, 'recording_list.mat');
info_file = fullfile(info_directory, 'session_info.mat');
kilosort_results_file = fullfile(kilosort_directory, 'kilosort_results.mat');

%%% Get recording list and session info
load(list_file, 'recording_list');
load(info_file, 'session_info');
n_recordings = session_info.n_recordings;
n_total_channels = session_info.n_total_channels;
start_frame_numbers = recording_list.start_frameno;
stop_frame_numbers = recording_list.stop_frameno;
stop_frame_numbers = stop_frame_numbers - start_frame_numbers(1);
start_frame_numbers = start_frame_numbers - start_frame_numbers(1);
sampling_rates = recording_list.sampling_rate;

%%% Get total length of all files if concatenated to potentially fill in
%%% gaps?
length_samples = stop_frame_numbers(end) + 1;
length_ms = ceil(length_samples * 1000) / sampling_rates(1);

%%% Get start times in ms
start_times = zeros(length(start_frame_numbers), 1);
for idx = 2:length(start_frame_numbers)
    start_times(idx) = start_frame_numbers(idx) * (1000 / sampling_rates(idx - 1));
end

%%% Get cumulative count of clusters and templates by recording, to give unique ID to
%%% each cluster. That way, if the session was divided into multiple
%%% recordings, and a channel or nearby channels were selected as the best
%%% channels for distinct clusters of distinct recordings, could then merge
%%% these clusters together based on similarities of the clusters.

cluster_count = 0; %%% Adding max value so that they don't repeat, given exclusion of some clusters, adding columns with new ID and original ID.
template_count = 0;

%%% Loop through recordings
for idx = 1:height(recording_list)
    
    %%% Get recording info
    recording_name = recording_list(idx, :).recording_name{:};
    n_channels = recording_list(idx, :).n_channels;
    n_samples = recording_list(idx, :).n_samples;
    sampling_rate = sampling_rates(idx);
    
    %%% Path with filtered recording to extract waveforms from
    recording_directory = fullfile(kilosort_directory, recording_name);
    load(fullfile(recording_directory, 'sorted_mapping.mat'), 'rowIDs');
    results_directory = fullfile(recording_directory, 'results');
    
    load(fullfile(results_directory, 'kilosort_results.mat'), 'kilosort_results');
    kilosort_configurations = kilosort_results.ops;
    nt0min = kilosort_configurations.nt0min; %%% Since kilosort padded the results, need to remove this amount of samples from start so that the waveforms are 2ms.
    
    %%% Load kilosort processed results for phy
    
    %%% Spike times in sample numbers
    spike_samples = readNPY(fullfile(results_directory, 'spike_times.npy'));
    %%% Convert spike samples to times in seconds
    spike_times = spike_samples/sampling_rate;
    
    %%% Spike amplitudes that were rescaled by kilosort
    spike_amplitudes = readNPY(fullfile(results_directory, 'amplitudes.npy'));
    
    %%% Template IDs for each spike
    spike_templates = readNPY(fullfile(results_directory, 'spike_templates.npy'));
    original_template_IDs = unique(spike_templates);
    new_template_IDs = original_template_IDs + template_count;
    template_count = template_count + length(original_template_IDs);
    
    %%% Cluster IDs fore each spike
    spike_clusters = readNPY(fullfile(results_directory, 'spike_clusters.npy'));
    original_cluster_IDs = unique(spike_clusters);
    new_cluster_IDs = original_cluster_IDs + cluster_count;
    cluster_count = cluster_count + length(original_cluster_IDs);
    
    
    channel_map = readNPY(fullfile(results_directory, 'channel_map.npy'));
    channel_map = channel_map + 1; %%% Convert from python to Matlab indices
    n_channels_map = length(channel_map);
    
    channel_positions = readNPY(fullfile(results_directory, 'channel_positions.npy'));
    
    templates = readNPY(fullfile(results_directory, 'templates.npy'));
    [n_templates, n_samples_waveforms, n_channels_templates] = size(templates);
    
    if any(spikes_templates ~= spikes_clusters)
        merge_map = get_merge_map(spikes_templates, spikes_clusters);
    else
        merge_map = [mat2cell(unique(original_cluster_IDs));mat2cell(unique(original_template_IDs))];
    end
    
    cluster_waveforms = 0 ; %%% For each cluster get the mean of the waveform of the templates that are part of the cluster
    best_template = 0;      %%% The one with the most counts for that cluster
    channel_ids = 0 ;       %%% Get the channels belonging to that template
    cluster_waveform;       %%% weighted based on counts of templates for each template
    
    spike_waveforms = 0 ;   %%% Extracting spike waveforms once best channels identified?
    
    cluster_amplitude = readtable(fullfile(results_directory, 'cluster_Amplitude.tsv'), 'FileType', 'text', 'Delimiter', '\t');
    cluster_contamination = readtable(fullfile(results_directory, 'cluster_ContamPct.tsv'), 'FileType', 'text', 'Delimiter', '\t');
    cluster_group = readtable(fullfile(results_directory, 'cluster_group.tsv'), 'FileType', 'text', 'Delimiter', '\t');
    cluster_kslabel = readtable(fullfile(results_directory, 'cluster_KSLabel.tsv'), 'FileType', 'text', 'Delimiter', '\t');
    similar_templates = readNPY(fullfile(results_directory, 'similar_templates.npy'));
    spike_clusters = readNPY(fullfile(results_directory, 'spike_clusters.npy'));
   
    whitening_mat = readNPY(fullfile(results_directory, 'whitening_mat.npy'));
    whitening_mat_inv = readNPY(fullfile(results_directory, 'whitening_mat_inv.npy'));
    a = 1;
end

end

function phy_processing()

%%% from .array import _index_of, _spikes_in_clusters, _spikes_per_cluster, SpikeSelector
%%% from .traces import (get_ephys_reader, RandomEphysReader, extract_waveforms, get_spike_waveforms, export_waveforms)

%%% def read_array(path, mmap_mode=None):
%%% Filter out nan and inf values.
%%% NOTE: do not check for nan/inf values on mmap arrays.
%%% TODO: virtual memmap array where the replacement is done on-the-fly when reading the array.

%%% def write_array(name, arr):
%%% Save an array to a binary file."""

%%% def from_sparse(data, cols, channel_ids):
%%% Convert a sparse structure into a dense one.
%%% Takes a matrix with n_spike rows and many columns representing different
%%% channels and only returns the data with columns representing the requested
%%% channel ids?

%%% def load_metadata(filename):
%%% Reads tsv amd returns rows with field name and columns being each cluster
%%% id and value

%%% def save_metadata(filename, field_name, metadata):
%%% Converts previous back to tsv

%%% def _all_positions_distinct(positions):
%%% Make sure all x y positions for channels are unique

%%% def get_closest_channels(channel_positions, channel_index, n=None):
%%% returns indices of n closest channels by euclidean distance, 1st being the channel index provided itself

%%% def _compute_pcs(x, npcs):
%%% returns n PCs of the spikes for each channel.

%%% def _project_pcs(x, pcs):
%%% projects the spikes into the pcs

%%% def compute_features(waveforms):
%%% takes care of above two

%%% def _find_first_existing_path(*paths, multiple_ok=True):
%%% There could be many paths, gets the first one that exists

%%% def _close_memmap(name, obj):
%%% ???


%%% Template model
% Spike attributes"
% -spike_clusters.npy
% -spike_templates.npy
% -spike_samples.npy
% -spike_times.npy
% -spike_times_reordered.npy
% -spike_amplitudes.npy

%%% template model object holds all data of kilosort dataset
templatemodel.dir_path = path_to_dataset;
templatemodel.dat_path = path_to_binary;
templatemodel.dtype = 'int16';
templatemodel.offset = 0;
templatemodel.n_channels = n_channels_in_recording;
templatemodel.sample_rate = sampling_rate;
n_closest_channels = 12; %%%Number of closest channels used for templates (average waveforms of clusters); hard coded in phy;
amplitude_threshold = 0; %%%Based on coding channels can be any amplitude relative to peak amplitude to be kept for templates
%%% Loads data
%%% Loads spike samples and spike times
%%% Makes sure the spike times are all in increasing order
       %%% "spike_times.npy" in samples or in seconds
%%% Loads amplitudes and checks that it is the same as number of spikes
       %%% 'amplitudes.npy', 'spikes.amps*.npy'
%%% Loads spike templates and makes sure first dim is number of spikes
       %%%'spike_templates.npy', 'spikes.templates*.npy'
       %%%Checks whether there are gaps in cluster IDs for templates
%%% Gets unique template_ids for spike_templates
%%% Loads spike_clusters and makes sure first dim is number of spikes
       %%%'spike_clusters.npy', 'spikes.clusters*.npy'
%%% Gets unique cluster_ids for spike clusters
%%% Gets spikes reordered if it exist and check that first dim is n spikes
%%% Gets channel_mapping and n_channels and makes sure that mapping dim is less
%%% than or equal to number of channels
       %%% 'channel_map.npy', 'channels.rawInd*.npy'
       %%%'channel_positions.npy', 'channels.localCoordinates*.npy'
%%% Loads channel shanks and makes sure dim is number of channels
       %%% 'channel_shanks.npy', 'channels.shanks*.npy'
%%% Loads channel probes and makes sure dim is number of channels
       %%%'channel_probe.npy', 'channels.probes*.npy'
%%% Loads sparse_templates makes sure dims  of data are n_templates, n_samples,
%%% n_channels, makes sure dims of cols is n_templates by n_channels
       %%%  'templates.npy', 'templates.waveforms.npy', 'templates.waveforms.*.npy'
       %%%# WARNING: KS2 saves templates_ind.npy (with an s), and note template_ind.npy,
       %%%# so that file is not taken into account here.
       %%%# That means templates.npy is considered as a dense array.
       %%%# Proper fix would be to save templates.npy as a true sparse array, with proper
       %%%# template_ind.npy (without an s).
       %%%path = self._find_path('template_ind.npy', 'templates.waveformsChannels*.npy')
       %%%cols = self._read_array(path)
       %%%if cols.ndim != 2:  # pragma: no cover
       %%%cols = np.atleast_2d(cols).T
       %%%assert cols.ndim == 2
       %%%logger.debug("Templates are sparse.")
       
       %%% data 'template_features.npy' cols 'template_feature_ind.npy' 'template_feature_spike_ids.npy'
       
%%% Loads spike_clusters and checks whether cluster number is different
%%% from template number to check the need to identify a merging map.
%%% Tries to load spike waveforms or otherwise fetches them from raw data
%%% Tries to identify whitening matrix
       %%% 'whitening_mat.npy' 'whitening_mat_inv.npy'
       %%% To unwhiten it is dot product of data and white m inv
%%% Loads similar_templates, template similiarity matrix n_templates x
%%% n_templates
       %%% 'similar_templates.npy'
%%% Loads traces and duration of recording and checks whether there are
%%% spike times that occur after end of recording
%%% Loads features
       %%% data 'pc_features.npy' cols 'pc_feature_ind.npy' rows 'pc_feature_spike_ids.npy'
%%% Loads sparse_features and determines n_features
%%%       Loads spike_attributes (any files that begin with spike_*.npy) turns
%%%       into dictionary
%%% Loads metadata
 

%%% Find best channels self, template, amplitude threshold
%%% Templates (n_samples, n_channels)
%%% Amplitude is determined as max-min of template for each channel
%%% Best channel is the channel with the largest amplitude
%%% Max amp amplitude of this channel
%%% Get indices of amplitudes that cross amplitude threshold * max_amp
%%% Find n closest channels, make sure that these channels cross threshold

%%% Template ID to get template to identify n best channels

%%% Might need to get template based on ID and unwhiten before getting best
%%% channels

%%% Get template sparse gets templates for a given id, removes channels
%%% with no signal by removing those that don't exceed 1e6 of the max only
%%% conserves templates and channels of a given id that did exceed that
%%% Also filters unused channels channel ids = -1. Unwhitens template if
%%% necessary. Tries to reorder channels_ids based on amplitude. Returns
%%% best channel, max amplitude, reordered template and channel ids based
%%% ona amplitude

%%% Map merges by getting template ids and getting cluster ids from their
%%% spikes and seeing how many clusters were in that one.

%%% Extracts waveforms from raw data

%%% Not necessary to compute PCs or project them?

%%% Depth based on spike pc features and weighted sum of coordinates
                %spikes_depths[ispi] = (np.sum(np.transpose(ypos * features) /
                 %                             np.sum(features, axis=1), axis=0))
%%% Sample2unit = 1 Uses templates to get amplitudes true, returns spike
%%% amplitudes, templates, average templates,  spike_amp = ks2_spike_amps * maxmin(inv_whitening(ks2_template_amps))
%%% np.matmul(data(n, :, :), wmi), amplitudes relative to min, templates_amps_au = np.max(templates_ch_amps, axis=1, spike_amps = templates_amps_au[spikes] * self.amplitudes
%%% template amplitude average of spikes of that template
%%% templates_amps_v = (np.bincount(spikes, weights=spike_amps) /np.bincount(spikes))
%%% templates_physical_unit = templates_wfs * (templates_amps_v / templates_amps_au)[:, np.newaxis, np.newaxis]
%%% return (spike_amps * sample2unit, templates_physical_unit * sample2unit, templates_amps_v * sample2unit)

%%% Get template from spikes. Spikes have template ids. Template ids are
%%% counted and sorted by number of spikes. Returns template with largest
%%% number of spikes from spike selection
end
