%%% MEA Pipeline Step 8: Compile MEA Recording Processing Information
%%%
%%% This function updates tables with processing information about MEA recording
%%% sessions and individual recordings by compiling the tables of each individual
%%% session in the database.

function n08_get_all_session_info(varargin)
if isempty(varargin)
    root_directory = '/path/to/MEA_pipeline/parent_directory';
else
    root_directory = varargin{1};  %%% (character array) Parent directory containing /MEA_pipeline and /MEA_database.
end

%%% Declare directories
data_directory = fullfile(root_directory, 'MEA_database');
session_list_file = fullfile(data_directory, 'all_sessions.mat');
recording_list_file = fullfile(data_directory, 'all_recordings.mat');

%%% Initialize output variables
all_sessions = [];
all_recordings = [];

%%% Search for individual session info and recording list files in database
info_directories = dir(fullfile(data_directory, '*', '*', 'info', 'session_info.mat'));

%%% If files are found put them all together
if ~isempty(info_directories)
    
    info_directories = {info_directories.folder};
    for idx = 1:length(info_directories)
        
        this_session = fullfile(info_directories{idx}, 'session_info.mat');
        these_recordings = fullfile(info_directories{idx}, 'recording_list.mat');
        
        if isfile(this_session)
            load(this_session, 'session_info');
        else
            session_info = [];
        end
        
        if isfile(these_recordings)
            load(these_recordings, 'recording_list');
        else
            recording_list = [];
        end
        
        all_sessions = [all_sessions; session_info];
        all_recordings = [all_recordings; recording_list];
    
    end

end

if ~isempty(all_sessions)
    save(session_list_file, 'all_sessions');
    save(recording_list_file, 'all_recordings');
end

end