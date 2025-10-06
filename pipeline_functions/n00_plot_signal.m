%%% The following function is for plotting raw and filtered recordings
%%% sample signal (1D and 2D views) and power spectral density
%%%
%%% Inputs:
%%% recordings: (n_channels x sampling_rate) matrix of amplitudes in uV
%%% root_directory: Directory of MEA project
%%% subjet: Subject ID
%%% folder: yyyy-mm-dd_name-of-session
%%% recording_name: eg. 'data0000', folder where recording stored
%%% signal_type: Whether signal is 'raw' or 'filtered'
%%% sampling_rate: recording sampling rate in Hz
%%%
%%% Output:
%%% For input recordings will make a single plot containing 3 subplots:
%%% -1D view of 1 second sample of signal of every channel (each in
%%% distinct color)
%%% -2D view of 1 second sample of signal colormapped by amplitude
%%% -Similar to first subplot, but for power spectral densities
%%% The figures are save in .png and .svg format in subject directory with
%%% name being concatenation of folder, recording name, and signal type.

function n00_plot_signal(recordings, root_directory, subject, folder, recording_name, signal_type, sampling_rate)

%%% Butterworth filter for extraction in combinato
min_pass = 300; %Hz
max_pass = 3000; %Hz

order = 6;
[butterworth_b, butterworth_a] = butter(order, [min_pass, max_pass] / (sampling_rate / 2), 'bandpass');

%%% Only keep a fraction of recorded channels
recordings = recordings(1:20:size(recordings, 1), :);

%%% Filter signal
recordings = filtfilt(butterworth_b, butterworth_a, double(recordings)')';

%%% Resort based on maximum absolute amplitude to highlight spikes
[~, sorted_indices] = sortrows(max(abs(recordings), [], 2), 'descend');
recordings = recordings(sorted_indices, :);

%%% Get recording matrix size
n_channels = size(recordings, 1);
n_samples = size(recordings, 2);

%%% Parameters for generating power spectra
frequency_resolution = sampling_rate / n_samples;           %%% Actual frequency resolution
min_frequency = 3;                                          %%% Minimum frequency for power spectrum plot
max_frequency = 606;                                        %%% Maximum frequency for power spectrum plot
frequencies = 0:frequency_resolution:(sampling_rate/2);     %%% The frequencies corresponding to fft result, 0 to nyquist freq
min_idx = find(frequencies >= min_frequency, 1, 'first');   %%% Indices to use to only store fft results of interest
max_idx = find(frequencies <= max_frequency, 1, 'last');

%%% Generate PSD and convert power to decibels
PSDs_recordings = 10 * log10(generate_PSD(recordings, sampling_rate, min_idx, max_idx));

%%% Initialize directory to hold plots for this subject
plot_directory = fullfile(root_directory, 'MEA_database', subject, 'plots');
if ~isfolder(plot_directory)
    mkdir(plot_directory);
end

plot_file = fullfile(plot_directory, strcat(folder, '_', recording_name, '_', signal_type));

%%% Plot parameters
figure_width = 1920;
figure_height = 1080; %%%Fit to screen 

subplot_width1 = round(figure_width * (3 / 5));
subplot_width2 = figure_width - subplot_width1;

subplot_height1 = round(figure_height / 7);
subplot_height2 = figure_height;
subplot_height3 = figure_height - subplot_height1;

xlim_signal = [1, n_samples];
xticks_signal = 0:(sampling_rate / 10):n_samples;
xticklabels_signal = arrayfun(@(x) sprintf('%.1fs', x), xticks_signal / sampling_rate, 'UniformOutput', false);
xticklabels_signal{1} = '';
xticklabels_signal{end} = '';

%%% Calculate appropriate limits to visualize signal. Tick marks spread from 0 by 50 uV.
max_absolute_amplitude = max(abs(recordings(:)));

if max_absolute_amplitude < 50

    ylim_signal = 50;
    yticks_signal = [-50, 0, 50];
    yticklabels_signal = {'-50uV', '0', '50'};
    
else

    ylim_signal = max_absolute_amplitude;
    yticks_signal = sort(unique([(-1 * ylim_signal), (ylim_signal), 0, 0:50:ylim_signal, 0:-50:(-1 * ylim_signal)]));
    yticklabels_signal = repelem({''}, 1, length(yticks_signal));
    yticklabels_signal{yticks_signal == -50} = '-50';
    yticklabels_signal{yticks_signal == 0} = '0';
    yticklabels_signal{yticks_signal == 50} = '50';
    
    if mod(max_absolute_amplitude, 50) >= 30 || mod(max_absolute_amplitude, 50) == 0
        yticklabels_signal{1} = strcat(num2str(round(-1 * max_absolute_amplitude)), 'uV');
        yticklabels_signal{end} = num2str(round(max_absolute_amplitude));
    else        
        yticklabels_signal{2} = strcat(num2str(round(yticks_signal(2))), 'uV');
        yticklabels_signal{end - 1} = num2str(round(yticks_signal(end - 1)));
    end
    
end
ylim_signal = [(-1.15 * ylim_signal), (1.15 * ylim_signal)]; %%% To make some extra space for xlabels

%%% Determine appropriate x and y ticks for PSD plots
xlim_PSD = [min_frequency, max_frequency];               %%% PSDs plotted in decibels (y) and linear frequency (x)
xticks_PSD = [min_frequency, 60:60:max_frequency];       %%% Each frequency tick centered at 60 Hz and harmonics
xticklabels_PSD = repelem({''}, 1, length(xticks_PSD));
xticklabels_PSD(2:2:end) = arrayfun(@(x) sprintf('%d', x), xticks_PSD(2:2:end), 'UniformOutput', false);
xticklabels_PSD{2} = '60Hz';
xticklabels_PSD{end} = '';

min_PSD = min(min(PSDs_recordings(~isinf(PSDs_recordings))));
max_PSD = max(max(PSDs_recordings(~isinf(PSDs_recordings))));

if max_PSD < 10 || isempty(max_PSD)
    max_PSD = 10;
end
if min_PSD > 0 || isempty(min_PSD)
    min_PSD = 0;
end  

ylim_PSD = [min_PSD, max_PSD];
yticks_PSD = sort(unique([min_PSD, max_PSD, 0, 0:10:max_PSD, 0:-10:min_PSD]));
ylim_PSD = [min(ylim_PSD)-(diff(ylim_PSD)*0.1), max(ylim_PSD)+(diff(ylim_PSD)*0.1)];
yticklabels_PSD = repelem({''}, 1, length(yticks_PSD));
yticklabels_PSD{yticks_PSD == 10} = '10dB';
yticklabels_PSD{yticks_PSD == 0} = '0';
if round(max_PSD) > 10
    yticklabels_PSD{yticks_PSD == max_PSD} = strcat(num2str(round(max_PSD)), 'dB');
end
if round(min_PSD) < 10 && round(min_PSD) ~= 0
    yticklabels_PSD{yticks_PSD == min_PSD} = strcat(num2str(round(min_PSD)), 'dB');
end

%%% Colormaps (first: would represent increase in row_ID, going from left to right columns in the MEA. 
%%% Used for signal and PSD. second: to represent amplitude of signal in 2D image.)
half_n = ceil(n_channels/2);

color_map_p1 = flipud(hot(ceil(half_n * 1.6)));
color_map_p2 = abyss(ceil(half_n * 1.6));
color_map = [color_map_p1(ceil(half_n * .6) + 1:end, :); color_map_p2(1:end - ceil(half_n * .6), :)];

%%% To highlight negative spikes with red to white color, and highlight
%%% positive spikes with blue, values closer to 0uV in black

%%% Y ticks and colorbar ticks for 2D signal view
ylim_2D = [0.5, n_channels + 0.5];
yticks_2D = 0:10:n_channels;
yticklabels_2D = [{''}, arrayfun(@(x) num2str(x), yticks_2D(2:end), 'UniformOutput', false)];
cticks_2D = yticks_signal;
cticklabels_2D = yticklabels_signal;

%%% Initialize figure to be the size of the screen and initialize matrix to
%%% hold axes for 3 subplots, the raw signal in 1D view, the raw signal in 2D view, and PSDs.

axes_handles = zeros(3, 1);

figure('Units', 'pixels', 'Visible', 'off')

figure_handle = gcf;
figure_handle.Position(3) = figure_width;
figure_handle.Position(4) = figure_height;

%%% Plot of raw signal, 1 second sample, slight offset and change in color
%%% for each individual channel
axes_handles(1, 1) = axes('Parent', figure_handle, 'Units', 'pixels', 'Position', [0, figure_handle.Position(4) - subplot_height1, subplot_width1, subplot_height1]);

hold on

for idx = 1:n_channels
    plot(recordings(idx, :), 'Color', color_map(idx, :), 'LineWidth', 1.25)
end

xlim(xlim_signal)
ylim(ylim_signal)
xticks(xticks_signal)
yticks(yticks_signal)
xticklabels([])
yticklabels([])

for idx = 1:length(xticks_signal)
    text(xticks_signal(idx), min(ylim_signal) + (0.08 * diff(ylim_signal)), xticklabels_signal{idx}, 'Color', [0 0 0], 'FontSize', 14, 'HorizontalAlignment', 'center');
end

for idx = 1:length(yticks_signal)

    if yticks_signal(idx) == 0
        text_color = [1 1 1];
    else
        text_color = [0 0 0];
    end
    
    text(sampling_rate * 0.015, yticks_signal(idx), yticklabels_signal{idx}, 'Color', text_color, 'FontSize', 14);
    
end

%%% Plot for PSDs
axes_handles(2, 1) = axes('Parent', figure_handle, 'Units', 'pixels', 'Position', [subplot_width1 + 1, 0, subplot_width2, subplot_height2]);

hold on

for idx = 1:n_channels
    plot(frequencies(min_idx:max_idx), PSDs_recordings(idx, :), 'Color', color_map(idx, :))
end

xlim(xlim_PSD)
ylim(ylim_PSD)
xticks(xticks_PSD)
yticks(yticks_PSD)
xticklabels([])
yticklabels([])

for idx = 1:length(xticks_PSD)
    if ~isempty(xticklabels_PSD{idx})
        text(xticks_PSD(idx), min(ylim_PSD)+(0.03 * diff(ylim_PSD)), xticklabels_PSD{idx}, 'Color', [0 0 0], 'FontSize', 16, 'HorizontalAlignment', 'center');
    end
end

for idx = 1:length(yticks_PSD)
    if ~isempty(yticklabels_PSD{idx})
        text(min(xlim_PSD) + (0.03 * diff(xlim_PSD)), yticks_PSD(idx), yticklabels_PSD{idx}, 'Color', [0 0 0], 'FontSize', 16);
    end
end

hold off

%%% Plot for 2D view of raw recording signal
axes_handles(3, 1) = axes('Parent', figure_handle, 'Units', 'pixels', 'Position', [0, 0, subplot_width1, subplot_height3]);

hold on

imagesc(1:n_samples, 1:n_channels, flipud(recordings));

colormap(color_map)

clim([-1 * max_absolute_amplitude, max_absolute_amplitude])

xlim(xlim_signal)
ylim(ylim_2D)
xticks(xticks_signal)
yticks(yticks_2D)
xticklabels([])
yticklabels([])

for idx = 1:length(xticks_signal)
    if ~isempty(xticklabels_signal{idx})
        text(xticks_signal(idx), min(ylim_2D) + (0.03 * diff(ylim_2D)), xticklabels_signal{idx}, 'Color', [1 1 1], 'FontSize', 16, 'HorizontalAlignment', 'center');
    end
end

for idx = 1:length(yticks_2D)
    if ~isempty(yticklabels_2D{idx})
        text(sampling_rate * 0.015, yticks_2D(idx), yticklabels_2D{idx}, 'Color', [1 1 1], 'FontSize', 16);
    end
end

colorbar_handle = colorbar('Location', 'southoutside');
colorbar_handle.YTick = cticks_2D;
colorbar_handle.YTickLabel = cticklabels_2D;
set(colorbar_handle, 'FontSize', 15)

hold off

%%% Print in .png for view and .svg for editing
print(plot_file, '-dpng');
print(plot_file, '-dsvg');
close all

end

%%% The following function generates power spectra for every channel (row) in
%%% data_samples and trims them to specified indices representing
%%% frequencies of interest 
function PSD = generate_PSD(data_samples, sampling_rate, min_idx, max_idx)
PSD = zeros(size(data_samples, 1), max_idx-min_idx+1);
N = size(data_samples, 2);
for idx = 1:size(data_samples, 1)
    sample = data_samples(idx, :);
    X = fft(sample); %Calculation of fourier coefficients
    Pxx = (abs(X(1:N/2+1)).^2) * 2 / (sampling_rate * N); %%% PSD calculation of real signal
    Pxx(1) = Pxx(1) / 2; %%% Correction of 0Hz DC value
    PSD(idx, :) = Pxx(min_idx:max_idx); %%% Only keeping the freq range to be plotted
end
end