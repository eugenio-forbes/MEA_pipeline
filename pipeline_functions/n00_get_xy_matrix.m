%%% Channel and electrode numbers of recording data
%%% from Maxwell Biosystems microelectrode arrays is not
%%% related to position of channel/electrode within the array.
%%% This function collects data from a microelectrode array to
%%% determine the pitch (distance between two electrodes) and
%%% based on the pitch it determines the number of channels in
%%% in the device. Then it produces vectors for all x and y
%%% coordinate pairs, arranged so that when reshaped this would
%%% result in matrices with the shape of the MEA. This would allow
%%% for remapping recording channels to row IDs so that recording
%%% data and results can more easily be manipulated and plotted.

function MEA_xy_matrix = n00_get_xy_matrix(varargin)
if isempty(varargin)
    root_directory = '/path/to/MEA_pipeline/parent_directory';
    subject = 'SC000';
    folder = 'yyyy-mm-dd_network';
else
    root_directory = varargin{1};  %%% (character array) Parent directory containing /MEA_pipeline and /MEA_database.
    subject = varargin{2};         %%% (character array) Subject code in 'SC000' format.
    folder = varargin{3};          %%% (character array) Folder name with data in 'yyyy-mm-dd_scan-name#' format. Raw recording in /folder/raw/data.raw.h5
end

%%% Declare directories
resource_directory = fullfile(root_directory, 'MEA_pipeline/base_code/resources');
data_directory = fullfile(root_directory, 'MEA_database');

%%% Get list of MEA devices, containing device ID in 1st column and pitch (distance of separation of electrodes in micrometers) in 2nd column.
device_list = fullfile(resource_directory, 'MEA_devices.mat');
if isfile(device_list)
    load(device_list, 'MEA_devices');
else
    MEA_devices = table;
end
    
%%% Get MEA device id corresponding to this session's raw file
maxwell_file = fullfile(data_directory, subject, folder, 'raw/data.raw.h5');
maxwell_struct = n00_maxwell_file_to_struct(maxwell_file);
MEA_id = maxwell_struct.wellplate.id{:};

if ~isempty(MEA_devices) && ismember(MEA_id, MEA_devices.id)
    %%% Get pitch for already listed device
    MEA_pitch = MEA_devices(ismember(MEA_devices.id, MEA_id), :).pitch;    
else
    %%% Determine pitch from mapping data and save new device id and
    %%% corresponding pitch to list
    recording_names = fields(maxwell_struct.data_store);
    n_recordings = length(recording_names);
    
    %%% Preallocate sufficiently large arrays of NaNs to hold all of the x
    %%% and y positions of the recordings.
    all_x = NaN(1024 * n_recordings, 1); %%% MEA only records from 1024 channels at a time
    all_y = NaN(1024 * n_recordings, 1);
    position_count = 0;
    
    for idx = 1:n_recordings
        this_x = maxwell_struct.data_store.(recording_names{idx}).settings.mapping.x;
        this_y = maxwell_struct.data_store.(recording_names{idx}).settings.mapping.yp;
        n_positions = length(this_x);
        all_x(position_count + 1:position_count + n_positions) = this_x;
        all_y(position_count + 1:position_count + n_positions) = this_y;
        position_count = position_count + n_positions;
    end
    
    pitch_x = min(diff(unique(all_x)));
    pitch_y = min(diff(unique(all_y)));
    MEA_pitch = min(pitch_x, pitch_y);
    
    if ~ismember(MEA_pitch, [17.5, 35])
        if MEA_pitch > 35
            error('Device %s didn''t gather enough data in %s %s recording to determine pitch.\n', MEA_id, subject, folder);
        elseif MEA_pitch < 17.5
            error('Device %s used in %s %s recording is likely a device with a higher resolution and code would need to be updated.\n', MEA_id, subject, folder);
        end
    else
        new_entry = table;
        new_entry.id = {MEA_id};
        new_entry.pitch = MEA_pitch;
        MEA_devices = [MEA_devices; new_entry];
        save(device_list, 'MEA_devices');
    end
end

%%% Based on pitch, set parameters for building MEA_xy_matrix
switch MEA_pitch
    case 17.5
        min_x = 0;
        max_x = 3832.5; 
        n_channels_x = 220;
        min_y = 0; 
        max_y = 2082.5; 
        n_channels_y = 120;
        total_channels = 26400;
    case 35
        min_x = 17.5; 
        max_x = 3832.5;
        n_channels_x = 110;
        min_y = 17.5;
        max_y = 2082.5;
        n_channels_y = 60;
        total_channels = 6600;
end

%%% Initialize MEA x and y matrix
MEA_xy_matrix = zeros(total_channels, 2); 
%%% First column x position value. Second column y position value.

unique_x = min_x:MEA_pitch:max_x;
unique_y = min_y:MEA_pitch:max_y;

all_x = repmat(unique_x, n_channels_y, 1);
all_y = repmat(unique_y', 1, n_channels_x);

MEA_xy_matrix(:, 1) = reshape(all_x, total_channels, 1);
MEA_xy_matrix(:, 2) = reshape(all_y, total_channels, 1);

%%% This way channels can be identified by row number instead of Maxlab's
%%% channel and electrode numbering system, so that any channel results can
%%% be input into a column vector, so that it can be reshaped/input into a
%%% matrix with the same dimensions as MEA and for localizing in
%%% plots without the need of x and y positions.

end