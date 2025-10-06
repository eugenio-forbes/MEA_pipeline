%%% MEA Pipeline, Step 7: Clean Up
%%%
%%% In development. Function meant for deletion of unnecessary files
%%% Once all other pipeline steps have been completed, such as int16
%%% copy of recording data.

function n07_clean_up(root_directory, subject, folder)
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
data_directory = fullfile(root_directory, 'MEA_database', subject, folder);
kilosort_directory = fullfile(data_directory, 'kilosort');

recording_directories = dir(kilosort_directory);
recording_directories = recording_directories(recording_directories.isdir);
recording_directories = fullfile({recording_directories.folder}, {recording_directories.name});

n_recordings = length(recording_directories;

for idx = 1:n_recordings
    recording_directory = recording_directories(idx)
    raw_recording_file = fullfile(recording_directory, 'raw_recording.bin');
    if isfile(raw_recording_file)
        delete(raw_recording_file)
    end
end

end