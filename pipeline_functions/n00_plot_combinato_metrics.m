%%% Function to plot metrics of Combinato results at the
%%% channel, cluster group, and cluster class levels

function n00_plot_combinato_metrics(varargin)
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

%%% Declare directories
data_directory = fullfile(root_directory, 'MEA_database', subject, folder);
info_directory = fullfile(data_directory, 'info');
kilosort_directory = fullfile(data_directory, 'kilosort');
spikes_directory = fullfile(data_directory, 'combinato_neurons');
if is_maxwell_data
    spikes_directory = [spikes_directory '_maxwell_data'];
end
matrices_directory = strrep(spikes_directory, 'combinato_neurons', 'matrices');

%%% Declare files
session_info_file = fullfile(info_directory, 'session_info.mat');
recording_list_file = fullfile(info_directory, 'recording_list.mat');
channels_file = fullfile(spikes_directory, 'channels.mat');
neurons_file = fullfile(spikes_directory, 'neurons.mat');
classes_file = fullfile(spikes_directory, 'classes.mat');

%%% Get info about MEA device used for recordings (pitch, number of channels, channels per axes), and recording names
load(session_info_file, 'session_info');
MEA_pitch = session_info.MEA_pitch;
n_total_channels = session_info.n_total_channels;
n_channels_x = session_info.n_channels_x;
n_channels_y = session_info.n_channels_y;
load(recording_list_file, 'recording_list');
recording_names = recording_list.recording_name;

%%% Load channels, neurons, classes
load(channels_file, 'channels');
load(neurons_file, 'neurons');
load(classes_file, 'classes');

neurons = struct2table(neurons);
classes = struct2table(classes);

%%% Declare excluded fields for channels
excluded_fields = {'row_ID', 'is_negative', 'subject', 'folder', 'n_recordings_detected', 'n_recordings_total', 'n_samples_total', 'duration_seconds_total'};

%%% Get channel table fields
channel_fields = channels.Properties.VariableNames;
channel_fields = channel_fields(~ismember(channel_fields, excluded_fields));

%%% Plot matrices for all fields in a single figure
%%% Separate by positive and negative spikes if not Maxwell data spikes (Maxwell only extracts negative spikes)
if is_maxwell_data
    is_negative = true;
    plot_channel_matrix(matrices_directory, channel_fields, channels, n_total_channels, n_channels_y, n_channels_x, MEA_pitch, is_negative);
else
    is_negative = false;
    positive_channels = channels(~channels.is_negative, :);
    plot_channel_matrix(matrices_directory, channel_fields, positive_channels, n_total_channels, n_channels_y, n_channels_x, MEA_pitch, is_negative);
    
    is_negative = true;
    negative_channels = channels(channels.is_negative, :);
    plot_channel_matrix(matrices_directory, channel_fields, negative_channels, n_total_channels, n_channels_y, n_channels_x, MEA_pitch, is_negative);
end

%%% Declare excluded fields for neurons
excluded_fields = {'row_ID', 'is_negative', 'subject', 'folder', 'recordings_detected', 'recordings_total', 'recordings_count', 'recordings_duration_seconds', 'group'};

%%% Get neuron table fields
neuron_fields = neurons.Properties.VariableNames;
neuron_fields = neuron_fields(~ismember(neuron_fields, excluded_fields));

%%% Plot histograms for all fields in a single figure
%%% Separate by positive and negative spikes if not Maxwell data spikes
plot_type = 'neurons';
if is_maxwell_data
    is_negative = true;
    plot_metric_histograms(matrices_directory, plot_type, neuron_fields, neurons, is_negative);
else
    is_negative = false;
    positive_neurons = neurons(~neurons.is_negative, :);
    plot_metric_histograms(matrices_directory, plot_type, neuron_fields, positive_neurons, is_negative);
    
    is_negative = true;
    negative_neurons = neurons(neurons.is_negative, :);
    plot_metric_histograms(matrices_directory, plot_type, neuron_fields, negative_neurons, is_negative);
end

%%% Declare excluded fields for classes
excluded_fields = {'row_ID', 'is_negative', 'subject', 'folder', 'recordings_detected', 'recordings_total', 'recordings_count', 'recordings_duration_seconds', 'group', 'class'};

%%% Get class table fields
class_fields = classes.Properties.VariableNames;
class_fields = class_fields(~ismember(class_fields, excluded_fields));

%%% Plot histograms for all fields in a single figure
%%% Separate by positive and negative spikes if not Maxwell data spikes
plot_type = 'classes';
if is_maxwell_data
    is_negative = true;
    plot_metric_histograms(matrices_directory, plot_type, class_fields, classes, is_negative);
else
    is_negative = false;
    positive_classes = classes(~classes.is_negative, :);
    plot_metric_histograms(matrices_directory, plot_type, class_fields, positive_classes, is_negative);
    
    is_negative = true;
    negative_classes = classes(classes.is_negative, :);
    plot_metric_histograms(matrices_directory, plot_type, class_fields, negative_classes, is_negative);
end

%%% Plot the channel maps for each individual recording in a single plot
plot_channel_maps(matrices_directory, kilosort_directory, recording_names, n_total_channels, MEA_pitch, n_channels_y, n_channels_x);

end

%%% For every metric of Combinato results, this function generates
%%% a color mapped surface contour displaying the results of every
%%% channel in the microelectrode array in the same shape as the MEA.
function plot_channel_matrix(matrices_directory, data_fields, channels, n_total_channels, n_channels_y, n_channels_x, MEA_pitch, is_negative)
n_fields = length(data_fields);
row_IDs = channels.row_ID;

if is_negative
    plot_name = 'MEA_metrics_negative_spikes';
else
    plot_name = 'MEA_metrics_positive_spikes';
end

%%% Plot parameters
figure_width = 1920;
figure_height = 1080;
title_space = 0;
xtick_space = 15;
ytick_space = 20;

if MEA_pitch < 35
    plot_xticks = 20:40:n_channels_x;
    plot_yticks = 20:40:n_channels_y;
else
    plot_xticks = 10:20:n_channels_x;
    plot_yticks = 10:20:n_channels_y;
end

plot_xticklabels = arrayfun(@(x) sprintf('%d', x * MEA_pitch), plot_xticks, 'UniformOutput', false);
plot_yticklabels = arrayfun(@(x) sprintf('%d', x * MEA_pitch), plot_yticks, 'UniformOutput', false);

%%% Calculate appropriate number of rows and columns and sizes, based on number of fields and sizes of margins to fit labels/titles.
n_columns = ceil(sqrt(n_fields));
n_rows = ceil(n_fields / n_columns);
axes_width = floor(figure_width / n_columns) - ytick_space;
axes_height = floor(figure_height / n_rows) - xtick_space - title_space;

%%% Initialize handles for axes
axes_handles = zeros(n_fields, 1);

figure('Units', 'pixels', 'Position', [0, 0, figure_width, figure_height], 'Visible', 'off');

for idx = 1:n_rows

    for jdx = 1:n_columns
    
        field_index = ((idx - 1) * n_columns) + jdx;
        
        if  field_index <= n_fields
            
            left_position = ((jdx - 1) * axes_width) + (jdx * ytick_space);
            top_position = figure_height - (idx * (axes_height + xtick_space + title_space));
            
            this_field = data_fields{field_index};
            
            plot_title = strrep(this_field, 'm_', 'mean ');
            plot_title = strrep(plot_title, 'std_', 'standard deviation ');
            plot_title = strrep(plot_title, '_', ' ');
            
            data_vector = zeros(n_total_channels, 1);
            data_vector(row_IDs) = channels{:, this_field};
            
            MEA_matrix = reshape(data_vector, n_channels_y, n_channels_x);
            
            axes_handles(field_index) = axes('Units', 'pixels', 'Position', [left_position, top_position, axes_width, axes_height]);
            
            imagesc(MEA_matrix);
            
            xticks(plot_xticks);
            yticks(plot_yticks);
            xticklabels(plot_xticklabels);
            yticklabels(plot_yticklabels);
            title(plot_title);
            
            ax = gca;
            ax.YAxis.TickLabelRotation = 90;
            set(gca, 'YDir', 'normal');
            
            colormap hot
            caxis([min(data_vector), max(data_vector)]);
            hcb = colorbar('Location', 'southoutside');
            set(hcb, 'FontSize', 12)
        end
        
    end
    
end

print(fullfile(matrices_directory, plot_name), '-dpng');
print(fullfile(matrices_directory, plot_name), '-dsvg');
close all

end

%%% Function to plot histograms for every metric of Combinato results
function plot_metric_histograms(matrices_directory, plot_type, data_fields, data_table, is_negative)

n_fields = length(data_fields);

if is_negative
    plot_name = sprintf('%s_metrics_negative_spikes', plot_type);
else
    plot_name = sprintf('%s_metrics_positive_spikes', plot_type);
end

%%% Plot parameters
figure_width = 1920;
figure_height = 1080;
title_space = 18;
xtick_space = 15;
ytick_space = 22;

%%% Calculate appropriate number of rows and columns and sizes, based on number of fields and sizes of margins to fit labels/titles.
n_columns = ceil(sqrt(n_fields));
n_rows = ceil(n_fields / n_columns);
axes_width = floor(figure_width / n_columns) - ytick_space;
axes_height = floor(figure_height / n_rows) - xtick_space - title_space;

%%% Initialize handles for axes
axes_handles = zeros(n_fields, 1);

figure('Units', 'pixels', 'Position', [0, 0, figure_width, figure_height], 'Visible', 'off');

for idx = 1:n_rows

    for jdx = 1:n_columns
    
        field_index = ((idx - 1) * n_columns) + jdx;
        
        if  field_index <= n_fields
            left_position = ((jdx - 1) * axes_width) + (jdx * ytick_space);
            top_position = figure_height - (idx * (axes_height + xtick_space + title_space)) + xtick_space;
            
            this_field = data_fields{field_index};
            
            plot_title = strrep(this_field, 'm_', 'mean ');
            plot_title = strrep(plot_title, 'std_', 'standard deviation ');
            plot_title = strrep(plot_title, '_', ' ');
            
            data_vector = data_table{:, this_field};
            if iscell(data_vector)
                data_vector = [data_vector{:}];
            end
            data_vector = data_vector(~isnan(data_vector) & ~isinf(data_vector));
            
            
            axes_handles(field_index) = axes('Units', 'pixels', 'Position', [left_position, top_position, axes_width, axes_height]);
            
            if ~isempty(unique(data_vector))
                if islogical(data_vector)
                    bin_edges = [-0.5, 0.5, 1.5];
                    histogram(data_vector, bin_edges);
                else
                    if length(unique(data_vector)) == 1
                        bin_edges = [unique(data_vector) - 0.5, unique(data_vector), unique(data_vector) + 0.5];
                        histogram(data_vector, bin_edges);
                    else
                        bin_width = range(data_vector)/100;
                        histogram(data_vector, 'BinWidth', bin_width);  
                    end
                end
            end
            
            title(plot_title);
            ax = gca;
            ax.YAxis.TickLabelRotation = 90;
        end
        
    end
    
end

print(fullfile(matrices_directory, plot_name), '-dpng');
print(fullfile(matrices_directory, plot_name), '-dsvg');
close all

end

%%% Function that for every recording in a recording session,
%%% will plot a color mapped surface contour displaying
%%% channels that were active during the recording.
function plot_channel_maps(matrices_directory, kilosort_directory, recording_names, n_total_channels, MEA_pitch, n_channels_y, n_channels_x)

n_recordings = length(recording_names);
plot_name = 'recording_channel_maps';

%%% Plot parameters
figure_width = 1920;
figure_height = 1080;
title_space = 18;
xtick_space = 15;
ytick_space = 20;

if MEA_pitch < 35
    plot_xticks = 20:40:n_channels_x;
    plot_yticks = 20:40:n_channels_y;
else
    plot_xticks = 10:20:n_channels_x;
    plot_yticks = 10:20:n_channels_y;
end

plot_xticklabels = arrayfun(@(x) sprintf('%d', x*MEA_pitch), plot_xticks, 'UniformOutput', false);
plot_yticklabels = arrayfun(@(x) sprintf('%d', x*MEA_pitch), plot_yticks, 'UniformOutput', false);

color_map = [[106, 133, 138] ; [217, 249, 165]] / 255; %%% Bluish gray for background and bright green for recording electrodes

%%% Calculate appropriate number of rows and columns and sizes, based on number of fields and sizes of margins to fit labels/titles.
n_columns = ceil(sqrt(n_recordings));
n_rows = ceil(n_recordings/n_columns);
axes_width = floor(figure_width/n_columns)-ytick_space;
axes_height = floor(figure_height/n_rows)-xtick_space-title_space;

%%% Initialize handles for axes
axes_handles = zeros(n_recordings, 1);

figure('Units', 'pixels', 'Position', [0, 0, figure_width, figure_height], 'Visible', 'off');

for idx = 1:n_rows

    for jdx = 1:n_columns
    
        field_index = ((idx - 1) * n_columns) + jdx;
        
        if  field_index <= n_recordings
        
            left_position = ((jdx - 1) * axes_width) + (jdx * ytick_space);
            top_position = figure_height - (idx * (axes_height + xtick_space + title_space)) + xtick_space;
            
            this_recording = recording_names{field_index};
            
            plot_title = this_recording;
            
            load(fullfile(kilosort_directory, this_recording, 'sorted_mapping.mat'), 'row_IDs')
            
            data_vector = zeros(n_total_channels, 1);
            data_vector(row_IDs) = 1;
            
            MEA_matrix = reshape(data_vector, n_channels_y, n_channels_x);
            
            axes_handles(field_index) = axes('Units', 'pixels', 'Position', [left_position, top_position, axes_width, axes_height]);
            
            imagesc(MEA_matrix);
            
            xticks(plot_xticks);
            yticks(plot_yticks);
            xticklabels(plot_xticklabels);
            yticklabels(plot_yticklabels);
            title(plot_title);
            
            ax = gca;
            ax.YAxis.TickLabelRotation = 90;
            set(gca, 'YDir', 'normal');
            
            colormap(color_map)
        end
        
    end
    
end

print(fullfile(matrices_directory, plot_name), '-dpng');
print(fullfile(matrices_directory, plot_name), '-dsvg');
close all

end