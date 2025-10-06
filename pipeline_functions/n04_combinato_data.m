%%% MEA Pipeline, Step 4: Gather Data from Combinato Results
%%%
%%% Function makes tables with spike clustering results at the
%%% channel, cluster group, and cluster class levels.

function n04_combinato_data(varargin)
if isempty(varargin)
    root_directory = '/path/to/MEA_pipeline/parent_directory';
    subject = 'SC000';
    folder = 'yyyy-mm-dd_network';
    is_maxwell_data = false;
else
    root_directory = varargin{1};  %%% (character array) Parent directory containing /MEA_pipeline and /MEA_database.
    subject = varargin{2};         %%% (character array) Subject code in 'SC000' format.
    folder = varargin{3};          %%% (character array) Folder name with data in 'yyyy-mm-dd_scan-name#' format. Raw recording in /folder/raw/data.raw.h5
    is_maxwell_data = varargin{4}; %%% (true/false)      Whether spike waveforms were extracted based on MW spike times.
end

warning off

%%% Declare directories and files
data_directory = fullfile(root_directory, 'MEA_database', subject, folder);
info_directory = fullfile(data_directory, 'info');
combinato_directory = fullfile(data_directory, 'combinato');
spike_directory = fullfile(data_directory, 'combinato_neurons');
if is_maxwell_data
    combinato_directory = [combinato_directory '_mw'];
    spike_directory = [spike_directory '_mw'];
    polarity_signs = {'neg'};
else
    polarity_signs = [{'neg'}; {'pos'}];
end
if ~isfolder(spike_directory)
    mkdir(spike_directory);
    
end
matrices_directory = strrep(combinato_directory, 'combinato', 'matrices');

%%% Load files with info about recordings and sorting
list_file = fullfile(info_directory, 'recording_list.mat');
sorted_channel_file = fullfile(matrices_directory, 'sorted_channels.mat');

load(list_file, 'recording_list');
load(sorted_channel_file, 'all_counts', 'all_indices', 'all_row_IDs');

%%% Loop through each channel/option/sign combination and convert .h5 file with results of sorting to
%%% matlab variables from which to extract spike data
%%% The results for each channel are stored in two .h5 files:
%%% 
%%% -data_(channel).h5:
%%%   -Single file with negative and positive spikes
%%%   -3 datasets grouped by /neg or /pos depending on polarity of spikes:
%%%
%%%     -/spikes              :       60 x n_spikes   : (single) 60 sample waveforms (2ms)
%%%     -/times               : n_spikes x 1          : (double) the times of n detected spikes
%%%     -/artifacts           : n_spikes x 1          : (int8) logical values pointing to artifactual spikes.
%%% 
%%% -sort_cat.h5:
%%%   -Separate files for clusterings of each channel and negative and positive spikes
%%%   -Does not match size of above array, some neurons skipped. 10 datasets:
%%%   
%%%     -/types               :        2 x n_groups   : (int16) for each group (1st row), the type (2nd row) (-1 artifact, 0 unclustered, 1 multi-unit, 2 single unit)
%%%     -/types_orig          :        2 x n_groups   : (int16) same as above, before manual clustering (not used)
%%%     -/artifacts           :        2 x n_classes  : (int64) for each class (1st row), whether it is artifact (1) or not (0) (2nd row)
%%%     -/artifacts_prematch  :        2 x n_classes  : (int64) same as above, before manual clustering (not used)
%%%     -/groups              :        2 x n_classes  : (int16) for each class (1st row), which group (2nd row)
%%%     -/groups_orig         :        2 x n_classes  : (int16) same as above, before manual clustering (not used)
%%%     -/classes             : (n_spikes-skipped) x 1: (uint16) class number for each spike
%%%     -/distance            : (n_spikes-skipped) x 1: (single) template matching distances
%%%     -/index               : (n_spikes-skipped) x 1: (uint32) python index (0:length(times)-1) of the spike times corresponding to spike (some skipped)
%%%     -/matches             : (n_spikes-skipped) x 1: (int8) 0 (SPC), 1 (1st template matching), 2 (second template matching)

%%% Loop through combinations of used row_IDs and sign
[Ax, Bx] = ndgrid(1:numel(polarity_signs), 1:numel(all_row_IDs));

polarity_signs = polarity_signs(Ax(:));

channels = table;
channels.row_ID = all_row_IDs(Bx(:));
n_rows = height(channels);

channels.is_negative = strcmp(polarity_signs, {'neg'});
channels.subject = repelem({subject}, n_rows, 1);
channels.folder = repelem({folder}, n_rows, 1);

%%% Initialize zero arrays the size of combination table for holding on to channel data
int8_zeros = int8(zeros(n_rows, 1));
int32_zeros = int32(zeros(n_rows, 1));
double_zeros = zeros(n_rows, 1);

channels.n_recordings_detected = int8_zeros;
channels.n_recordings_total = int8_zeros;
channels.n_samples_total = int32_zeros;
channels.dur_seconds_total = double_zeros;
channels.n_spike_detections = int32_zeros;
channels.n_good_groups = int8_zeros;
channels.n_bad_groups = int8_zeros;
channels.n_good_classes = int16_zeros;
channels.n_bad_classes = int16_zeros;
channels.n_good_spikes = int32_zeros;
channels.n_bad_spikes = int32_zeros;
channels.n_skipped_spikes = int32_zeros;
channels.n_spc_match = int32_zeros;
channels.n_1st_template_match = int32_zeros;
channels.n_2nd_template_match = int32_zeros;
channels.n_sua = int8_zeros;
channels.n_mua = int8_zeros;
channels.n_noise = int8_zeros;
channels.isi_SNR = int8_zeros; 
channels.m_firing_rate = double_zeros;
channels.m_amplitude = double_zeros;

%%% Intitialize empty cell array to hold on to group/unit/neuron data.
neurons = cell(n_rows, 1);

%%% Intitialize empty cell array to hold on to subgroup/classes data.
classes = cell(n_rows, 1);

%%% One logical array to keep track of which combinations actually had
%%% spike data
no_data = true(n_rows, 1);

%%% Initialize parallel pool if it hasn't been created/if it time out.
poolobj = gcp('nocreate');
if ~isempty(poolobj)
    delete(poolobj);
end
parpool(36)

parfor idx = 1:n_rows
    
    row_ID = channels(idx, :).row_ID;
    is_negative = channels(idx, :).is_negative;
    
    this_folder = sprintf('MEA_%05d', row_ID);
    spike_data_file = fullfile(combinato_directory, this_folder, sprintf('data_MEA_%05d.h5', row_ID));
    sort_file = fullfile(combinato_directory, this_folder, sprintf('sort_%s_options00', polarity_signs{idx}), 'sort_cat.h5');
    
    if isfile(spike_data_file) && isfile(sort_file) %%% Because not every row_ID had positive or negative spikes detected
        
        no_data(idx) = false;
        
        %%% Get info specific to channel
        if is_negative
            if is_maxwell_data
                sign_idx = 1;
            else
                sign_idx = 2;
            end
        else
            sign_idx = 1;
        end
        
        channels(idx, :).n_recordings_detected = all_counts(row_ID, sign_idx);
        recording_indices = find(all_indices(row_ID, :));
        channels(idx, :).n_recordings_total = length(recording_indices);
        
        recording_samples = recording_list(recording_indices, :).n_samples;
        channels(idx, :).n_samples_total = sum(recording_samples);
        
        recording_sampling_rates = recording_list(recording_indices, :).sampling_rate;
        channels(idx, :).dur_seconds_total = sum(recording_samples./recording_sampling_rates);
        
        %%% Read variables from spike data file and sort file    
        spikes = struct('waveforms', [], 'times', [], 'artifacts', []);
        spikes.waveforms = h5read(spike_data_file, sprintf('/%s/spikes', polarity_signs{idx}));
        spikes.times = h5read(spike_data_file, sprintf('/%s/times', polarity_signs{idx}));
        spikes.artifacts = h5read(spike_data_file, sprintf('/%s/artifacts', polarity_signs{idx}));
        
        sorting_results = struct('types', [], 'artifacts', [], 'groups', [], 'classes', [], 'index', [], 'matches', [], 'distances', []);
        sorting_results.types = h5read(sort_file, '/types');
        sorting_results.artifacts = h5read(sort_file, '/artifacts');
        sorting_results.groups = h5read(sort_file, '/groups');
        sorting_results.classes = h5read(sort_file, '/classes');
        sorting_results.index = h5read(sort_file, '/index');
        sorting_results.matches = h5read(sort_file, '/matches');
        sorting_results.distances = h5read(sort_file, '/distance');
        
        %%% Make channel_info structure for use in building neurons and classes tables
        channel_info = struct;
        channel_info.subject = subject;
        channel_info.folder = folder;
        channel_info.row_ID = row_ID;
        channel_info.is_neg = is_negative;
        channel_info.recordings_detected = channels(idx, :).n_recordings_detected;
        channel_info.recordings_count = channels(idx, :).n_recordings_total;
        channel_info.recordings_dur_seconds = channels(idx, :).dur_seconds_total;
        
        %%% Channel info logged, artifacts removed, spikes ordered by group then index, with group number (1) start (2) and end (3) indices provided.
        [this_channel, ~, ISIs, waveforms, matches, distances, spike_classes, class_groups, group_indices] = transform_spike_data(channel_info, spikes, sorting_results);
        
        %%% Log data to channel table
        channels(idx, :).n_spike_detections = this_channel.n_detections;
        channels(idx, :).n_good_groups = this_channel.n_groups;
        channels(idx, :).n_bad_groups = this_channel.n_bad_groups;
        channels(idx, :).n_good_classes = this_channel.n_classes;
        channels(idx, :).n_bad_classes = this_channel.n_bad_classes;
        channels(idx, :).n_good_spikes = this_channel.n_spikes;
        channels(idx, :).n_bad_spikes = this_channel.n_bad_spikes;
        channels(idx, :).n_skipped_spikes = this_channel.n_skipped;
        channels(idx, :).n_spc_match = this_channel.n_spc;
        channels(idx, :).n_1st_template_match = this_channel.n_1st;
        channels(idx, :).n_2nd_template_match = this_channel.n_2nd;
        
        %%% Development:Could save artifact free and group-ordered variables in a lighter format, 
        %%% but for now will keep data in combinato format and only log information of each hierarchical
        %%% level for use in accessing and refining data in analysis.
        
        if this_channel.n_classes > 0 && this_channel.n_groups > 0
            %%% Classes structure created, class info logged, class silhouette
            %%% scores calculated. Classes structure created before neuron
            %%% structure to pass as an argument silhouette scores
            these_classes = make_class_structure(channel_info, this_channel, ISIs, waveforms, matches, distances, spike_classes, class_groups);
            
            %%% Neurons (group) structure created, group info logged
            [this_channel, these_neurons] = make_neuron_structure(channel_info, this_channel, ISIs, waveforms, matches, distances, spike_classes, class_groups, group_indices);
            
            %%% Log channel information about groups
            channels(idx, :).n_sua = this_channel.n_sua;
            channels(idx, :).n_mua = this_channel.n_mua;
            channels(idx, :).n_noise = this_channel.n_noise;
            channels(idx, :).isi_SNR = this_channel.isi_SNR;
            channels(idx, :).m_firing_rate = sum([these_neurons.mean_FR].*[these_neurons.n_spikes])/sum([these_neurons.n_spikes]);
            channels(idx, :).m_amplitude = sum([these_neurons.mean_peak_uV].*[these_neurons.n_spikes])/sum([these_neurons.n_spikes]);

            
            %%% Concatenating structures associated to a session's bank (sorted_micros)
            neurons{idx} = these_neurons;
            classes{idx} = these_classes;
        end
    end
end

%%% Close parallel pool
if ~isempty(poolobj)
    delete(poolobj);
end

neurons = vertcat(neurons{:});
classes = vertcat(classes{:});

save(fullfile(spike_directory, 'channels.mat'), 'channels');
save(fullfile(spike_directory, 'neurons.mat'), 'neurons');
save(fullfile(spike_directory, 'classes.mat'), 'classes');

end


function [channel, times, ISIs, waveforms, spike_matches, spike_distances, spike_classes, class_groups, group_indices] = transform_spike_data(channel_info, spikes, sorting_results)

old_times = spikes.times; 
n_detections = length(old_times);
old_waveforms = double(spikes.waveforms);
spike_artifacts = spikes.artifacts;
spike_artifacts = find(spike_artifacts);
group_types = double(sorting_results.types);
class_groups = double(sorting_results.groups);
class_artifacts = double(sorting_results.artifacts);
class_artifacts = class_artifacts(1, class_artifacts(2, :)==1);
spike_classes = double(sorting_results.classes);
spike_indices = double(sorting_results.index);
spike_matches = double(sorting_results.matches);
spike_distances = double(sorting_results.distances);

n_skipped = size(old_waveforms, 2) - length(spike_indices);
bad_types = group_types(2, :)<1;
bad_groups = group_types(1, bad_types);
bad_classes = class_groups(1, ismember(class_groups(2, :), bad_groups));
bad_classes = unique([bad_classes, class_artifacts]);
bad_indices = spike_indices(ismember(spike_classes, bad_classes));
bad_indices = unique([bad_indices;(spike_artifacts-1)]); 

%%% Combinato uses python indices, so subtract 1 from artifacts to find good indices

%%% Doing it this way because some spike indices are skipped
good_spikes = ~ismember(spike_indices, bad_indices);
good_indices = spike_indices(good_spikes);
spike_classes = spike_classes(good_spikes);
spike_matches = spike_matches(good_spikes);
spike_distances = spike_distances(good_spikes);
class_groups = class_groups(:, ismember(class_groups(1, :), unique(spike_classes)));
table_groups = arrayfun(@(x) class_groups(2, class_groups(1, :)==x), spike_classes);

%%% Make table to sort spike data by group then index
index_table = table;
index_table.index = good_indices;
index_table.class = spike_classes;
index_table.group = table_groups;
index_table.distance = spike_distances;
index_table.match = spike_matches;
index_table = sortrows(index_table, {'group', 'index'}, 'ascend');

spike_classes = index_table.class;
spike_distances = index_table.distance;
spike_matches = index_table.match;
table_groups = index_table.group; 
good_indices = index_table.index;
clear index_table;

good_indices = good_indices + 1; %%% Add + 1 to keep only the good spikes.

times = old_times(good_indices); clear old_times %%%Sorted by group and index
waveforms = old_waveforms(:, good_indices);
waveforms = waveforms';
ISIs = diff(times);
ISIs(ISIs<0)= NaN;

good_classes = unique(spike_classes, 'stable');
good_groups = unique(table_groups, 'stable');

group_indices = arrayfun(@(x) [x, find(table_groups==x, 1, 'first'), find(table_groups==x, 1, 'last')], good_groups, 'UniformOutput', false);
group_indices = cell2mat(group_indices);

%%% Create channel structure with info about clustering results
channel = channel_info;
channel.n_detections = n_detections;
channel.n_groups = length(good_groups);
channel.n_bad_groups = length(bad_groups);
channel.n_classes = length(good_classes);
channel.n_bad_classes = length(bad_classes);
channel.n_spikes = length(times);
channel.n_bad_spikes = length(bad_indices);
channel.n_skipped = n_skipped;
channel.n_features = size(waveforms, 2);
channel.n_spc = sum(spike_matches == 0);
channel.n_1st = sum(spike_matches == 1);
channel.n_2nd = sum(spike_matches == 2);

end


function classes = make_class_structure(channel_info, channel, ISIs, waveforms, matches, distances, spike_classes, class_groups)
%%% Get channel info
n_classes = channel.n_classes;
n_groups = channel.n_groups;
n_features = channel.n_features;
n_detections = channel.n_detections;
is_neg = channel.is_neg;

%%% Define line noise interspike intervals in ms, floored for indexing in histogram
harmonics_ISIs = floor(1000./(60.*(1:5))); %in ms the expected ISIs of 60Hz to 5th harmonic
pseudoSubH_ISIs = floor(1000./(60./(2:15))); %in ms ISIs of seeming subharmonics occuring in multiples of 16.66ms
line_noise_ISIs = [harmonics_ISIs, pseudoSubH_ISIs];

%%% Create array with individual spike group labels
unique_classes = class_groups(1, :)';

%%% Start classes structure from channel_info
classes = repmat(channel_info, n_classes, 1);

%%% Loop through classes
for kdx = 1:n_classes
    %%% Get class and group number
    this_class = unique_classes(kdx);
    this_group = class_groups(2, unique_classes==this_class);
    
    %%% Find indices of spikes in the class
    idx = find(spike_classes == this_class);
    n_spikes = length(idx); %n_spikes in this class
    
    %%% Get spike data
    class_matches = matches(idx);
    class_distances = distances(idx);
    
    %%% For ISI indices subtract one and remove 0 index (first spike not counted)
    idxISIs = idx-1;
    idxISIs = idxISIs(idxISIs~=0);
    class_ISIs = ISIs(idxISIs); %ISIs are as calculated within entire group.
 
    %%% Calculate percentage of ISIs below 3 ms and contributing to line noise
    p_sub3ms = mean(class_ISIs<3); %Proportion of spikes in this class that contribute to <3 ms in group
    hist_ISIs = histcounts(class_ISIs, 1:256); %Not counting 0 to 1 for indices to match edges
    p_60Hz = sum(hist_ISIs(line_noise_ISIs))/n_spikes; %Proportion of spikes in this class that contribute to line noise ISIs
        
    %%% Get class absolute waveforms for many other metrics
    class_waveforms = waveforms(idx, :);
    if is_neg
        class_waveforms = class_waveforms * -1; %In order to get absolute peak values
    end
    
    %%% Log class info
    
    %%% Class and group numbers
    classes(kdx).class = this_class;
    classes(kdx).group = this_group;
    
    %%% The total number of spikes in the class, and as a faction from detected in the channel recording
    classes(kdx).n_spikes = n_spikes;                           
    classes(kdx).f_spikes = n_spikes/n_detections;
    
    %%% The percentage of spikes in this class that contribute to sub3ms and line noise ISIs
    classes(kdx).p_sub3ms = p_sub3ms;
    classes(kdx).p_60Hz = p_60Hz;
    
    %%% Percentages in this class that were matched through SPC, 1st template or 2nd template matching
    classes(kdx).p_spc = sum(class_matches == 0)/length(class_matches);    
    classes(kdx).p_1st = sum(class_matches == 1)/length(class_matches);
    classes(kdx).p_2nd = sum(class_matches == 2)/length(class_matches);
    
    %%% Mean and standard deviation of the class's spike amplitudes in uV (absolute value).
    classes(kdx).m_amp_uV = mean(max(class_waveforms, [], 1));              
    classes(kdx).std_amp_uV = std(max(class_waveforms, [], 1));
    
    %%% Mean and standard deviation of average voltage of class's waveforms
    classes(kdx).m_voltage = mean(mean(class_waveforms, 1));                
    classes(kdx).std_voltage = std(mean(class_waveforms, 1));               
    
    %%% Mean and standard deviation of all matching distances
    classes(kdx).m_distances = mean(class_distances);                      
    classes(kdx).std_distances = std(class_distances);
end

end

function [channel, neurons] = make_neuron_structure(channel_info, channel, ISIs, waveforms, matches, distances, spike_classes, class_groups, group_indices)
%%% Get channel info
n_groups = channel.n_groups;
n_detections = channel.n_detections;
n_features = channel.n_features;
is_neg = channel.is_neg;
dur_seconds = channel.recordings_dur_seconds;

%%% Define line noise interspike intervals in ms, floored for indexing in histogram
harmonics_ISIs = floor(1000./(60.*(1:5))); %in ms the expected ISIs of 60Hz to 5th harmonic
pseudoSubH_ISIs = floor(1000./(60./(2:15))); %in ms ISIs of seeming subharmonics occuring in multiples of 16.66ms
line_noise_ISIs = [harmonics_ISIs, pseudoSubH_ISIs];
not_noise_ISIs = ~ismember(1:255, line_noise_ISIs);

%%% Start neurons (group) structure from channel_info
neurons = repmat(channel_info, n_groups, 1);

%%% Make a spike group label array
unique_classes = class_groups(1, :);
spike_groups = arrayfun(@(x) class_groups(2, unique_classes == x), spike_classes);
unique_groups = group_indices(:, 1);

%%% Loop through groups
for kdx = 1:n_groups
    
    %%% Get indices of spikes for the group
    this_group = unique_groups(kdx);
    idx = find(spike_groups == this_group);
    
    %%% Get number of spikes and firing rate
    n_spikes = length(idx);
    firing_rate = n_spikes/dur_seconds;
    
    %%% Determine how many classes in group
    group_classes = class_groups(2, :)==this_group;
    n_classes = sum(group_classes);
    
    %%% Get spike data
    group_spike_matches = matches(idx);
    group_spike_distances = distances(idx);
    % For ISIs subtract 1 from idx and remove 0 index (first spike)
    ISI_idx = idx-1;
    ISI_idx = ISI_idx(ISI_idx~=0);
    group_ISIs = ISIs(ISI_idx);
    
    %%% Calculate percentage of ISIs below 3 ms and contributing to line noise
    p_sub3ms = mean(group_ISIs<3); %Percentage of ISIs in group < 3ms
    hist_ISIs = histcounts(group_ISIs, 1:256); %Not counting 0 to 1 for indices to match edges
    p_60Hz = sum(hist_ISIs(line_noise_ISIs))/n_spikes; %Percentage of ISIs in line noise ISIs
    
    %%% Calculate ratio of mean line_noise counts to not noise counts
    isi_SNR = mean(hist_ISIs(not_noise_ISIs))/mean(hist_ISIs(line_noise_ISIs));
    
    group_spikes = waveforms(idx, :);
    if is_neg
        group_spikes = group_spikes * -1; %In order to get absolute peak values
    end
    centroid = mean(group_spikes, 1);
    group_std_spikes = std(group_spikes, [], 1);
    
    %%% Get mean waveform peak, width, peak std and avg std of last 25% samples
    mean_peak_uV = max(centroid);
    mean_peak_idx = find(centroid == mean_peak_uV, 1, 'first');
    mean_trough = min(centroid(1:mean_peak_idx-1));
    mean_trough_idx = find(centroid(1:mean_peak_idx-1) == mean_trough, 1, 'first');
    peak_width =  mean_peak_idx - mean_trough_idx + 1;
    std_peak = std(group_spikes(:, mean_peak_idx));
    std_last25 = mean(group_std_spikes(round(n_features*0.75):end));
    inflection_points = find(diff(diff(centroid))==0);
    inflection_points = inflection_points(inflection_points>mean_peak_idx-2);
    if ~isempty(inflection_points)
        inflection_values = [mean_peak_uV, centroid(inflection_points)];
        inflection_magnitudes = abs(diff(inflection_values));
        num_inflection_violations = sum(inflection_magnitudes > (0.2*abs(mean_peak_uV)));
    else
        num_inflection_violations = 0;
    end
    
    %%% Get group grade based on criteria
    grade = grade_group(firing_rate, p_sub3ms, p_60Hz, isi_SNR, mean_peak_uV, std_last25, std_peak, num_inflection_violations);
    
    %%% Log class info
    %%% Group number
    neurons(kdx).group = this_group;
    
    %%% Grade divided in 3 bins to filter out a certain type
    neurons(kdx).is_SUA = strcmp(grade, 'SUA');
    neurons(kdx).is_MUA = strcmp(grade, 'MUA');
    neurons(kdx).is_noise = strcmp(grade, 'noise');
    
    %%% The total number of spikes in the class, and as a faction from detected in the channel recording
    neurons(kdx).n_spikes = n_spikes;                           
    neurons(kdx).f_spikes = n_spikes/n_detections;
    neurons(kdx).n_classes = n_classes;
    
    %%% The mean firing rate in Hz
    neurons(kdx).mean_FR = firing_rate;
    neurons(kdx).mean_peak_uV = mean_peak_uV;
    neurons(kdx).peak_width = peak_width;
    
    %%% The percentage of spikes in this class that contribute to sub3ms and line noise ISIs
    neurons(kdx).p_sub3ms = p_sub3ms;
    neurons(kdx).p_60Hz = p_60Hz;
    neurons(kdx).isi_SNR = isi_SNR;
    
    %%% Percentages in this class that were matched through SPC, 1st template or 2nd template matching
    neurons(kdx).p_spc = sum(group_spike_matches == 0)/length(group_spike_matches);    
    neurons(kdx).p_1st = sum(group_spike_matches == 1)/length(group_spike_matches);
    neurons(kdx).p_2nd = sum(group_spike_matches == 2)/length(group_spike_matches);
    
    %%% Mean and standard deviation of the group's spike amplitudes in uV (absolute value).
    neurons(kdx).m_amp_uV = mean(max(group_spikes, [], 1));              
    neurons(kdx).std_amp_uV = std(max(group_spikes, [], 1));
    
    %%% Mean and standard deviation of average voltage of class's waveforms
    neurons(kdx).m_voltage = mean(mean(group_spikes, 1));                
    neurons(kdx).std_voltage = std(mean(group_spikes, 1));               
    
    %%% Mean and standard deviation of all matching distances
    neurons(kdx).m_distances = mean(group_spike_distances);                      
    neurons(kdx).std_distances = std(group_spike_distances);
end

%%% Get channel totals

channel.n_sua = sum([neurons.is_SUA]);
channel.n_mua = sum([neurons.is_MUA]);
channel.n_noise = sum([neurons.is_noise]);
channel.isi_SNR = mean([neurons.isi_SNR]);    

end


function grade = grade_group(firing_rate, p_sub3ms, p_60Hz, isi_SNR, mean_peak_uV, std_last25, std_peak, num_inflection_violations)
%%% Grading criteria adapted from doi:10.1088/1741-2560/10/1/016001
%%% Substituted PSD criteria for other line noise features

grade = 'noise'; %until proven otherwise

%%% Criterion 1: Firing rate above a threshold
fr_sua = firing_rate > 0.05;
fr_mua = firing_rate > 0.05;

%%% Criterion 2: Refractory period violation below a threshold
sub3ms_sua = p_sub3ms < 0.05;
sub3ms_mua = p_sub3ms < 0.1;

%%% Criterion 3: Ratio of peak amplitude to standard deviation of last 25% samples
peak_ratio = mean_peak_uV/std_last25;
peak_ratio_sua = peak_ratio > 1;
peak_ratio_mua = peak_ratio > 0.5;

%%% Criterion 4: Ratio of peak std to stds of last 25% samples
std_ratio = std_peak/std_last25;
std_ratio_sua = std_ratio>0.33;
std_ratio_mua = std_ratio>0.33;

%%% Criterion 5: Number of inflection magnitude violations
inflections_sua = num_inflection_violations <3;
inflections_mua = num_inflection_violations <= 3;

%%% Must meet all criteria to be graded mua or sua
if fr_mua && sub3ms_mua && peak_ratio_mua && std_ratio_mua && inflections_mua
    grade = 'MUA';
end
if fr_sua && sub3ms_sua && peak_ratio_sua && std_ratio_sua && inflections_sua
    grade = 'SUA';
end

%%% Adding this to criteria. If exceeds certain noise levels then grade as noise instead
if p_60Hz > 0.2 || isi_SNR < 0.8
    grade = 'noise';
end

end