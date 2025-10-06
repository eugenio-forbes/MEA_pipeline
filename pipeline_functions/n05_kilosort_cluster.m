%%% MEA Pipeline, Step 5: Using Kilosort to Cluster Spikes
%%%
%%% NOTE! This function requires CUDA and must be executed from
%%% a GPU node.
%%%
%%% This function executes signal processing, spike extraction,
%%% and spike clustering using Kilosort functions.
%%%
%%% NOTE! Two Kilosort functions, datashift() and rezToPhy()
%%% have been modified to be compatible with MEA data.

function n05_kilosort_cluster(varargin)
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
info_directory = fullfile(data_directory, 'info');
error_directory = fullfile(root_directory, 'MEA_database/error_logs', subject, folder);

%%% Get recording list and session info
list_file = fullfile(info_directory, 'recording_list.mat');
info_file = fullfile(info_directory, 'session_info.mat');
load(list_file, 'recording_list');
load(info_file, 'session_info');
n_recordings = session_info.n_recordings;
MEA_pitch = session_info.MEA_pitch;

for idx = 1:n_recordings
    %%% Gather info for each recording
    recording_info = recording_list(idx, :);
    recording_name = recording_info.recording_name{:};
    
    %%% To save kilosort results visualizable in Phy
    recording_directory = fullfile(kilosort_directory, recording_name);
    results_directory = fullfile(recording_directory, 'results');
    if ~isfolder(results_directory)
        mkdir(results_directory);
    end
    
    %%% Use this info to tailor congigurations to be used by kilosort
    %%% Other parameters can be modified within this function.
    kilosort_configurations = n00_set_kilosort_configs(recording_directory, recording_info, MEA_pitch);
    
    try
        
        %%% All the functions used by kilosort below
        %%% These functions have been modified to adapt to MEA, make plots
        %%% invisible, and save them to folder. Also disabled all displays
        %%% of messages except for those written in this function to track
        %%% progress.
        %%% Results of processing and sorting saved in rez and .npy files for Phy (results)
 
        %%% Data filtering and 'whitening'
        fprintf('%s : Starting processing of raw recording in %s for subject %s, %s session.\n', datestr(now), recording_name, subject, folder);
        kilosort_results = preprocessDataSub(kilosort_configurations);
        
        %%% Extraction of spikes and assignment of them to 'upsampled' x
        %%% and y positions. Mainly for getting upsampled positions???
        fprintf('%s : Upsampling x and y positions.\n', datestr(now));
        kilosort_results = datashift2(kilosort_results, 1, recording_directory, MEA_pitch);
        
        %%% Actual extraction of spikes with upsampled positions, using
        %%% closest channels to the upsampled positions
        fprintf('%s : Extracting spikes.\n', datestr(now));
        [kilosort_results, st3, tF] = extract_spikes(kilosort_results);
        
        %%% Template learning
        fprintf('%s : Performing template learning.\n', datestr(now));
        kilosort_results = template_learning(kilosort_results, tF, st3);
        
        %%% Sorting
        fprintf('%s : Sorting spikes.\n', datestr(now));
        [kilosort_results, st3, tF] = trackAndSort(kilosort_results);
        
        %%% Clustering
        fprintf('%s : Clustering and merging spikes.\n', datestr(now));
        kilosort_results = final_clustering(kilosort_results, tF, st3);
        
        %%% Merging
        fprintf('%s : Merging clusters.\n', datestr(now));
        kilosort_results = find_merges(kilosort_results, 1);
        
        %%% Saving results to .npy files used by Phy
        fprintf('Saving kilosort clustering results.\n');
        rezToPhy2(kilosort_results, results_directory);
        
        %%% Saving rez structure to results directory as well
        save(fullfile(results_directory, 'kilosort_results.mat'), 'kilosort_results');
        
        fprintf('%s : Finished kilosorting recording in %s for subject %s, %s session.\n', datestr(now), recording_name, subject, folder);

        
    catch this_error
        fprintf('%s : ERROR kilosorting recording in %s for subject %s, %s session.\n', datestr(now), recording_name, subject, folder);
        this_error_directory = fullfile(error_directory, recording_info.recording_name{:});
        if ~isfolder(this_error_directory)
            mkdir(this_error_directory);
        end
        error_file = fullfile(this_error_directory, 'n05_kilosort_cluster_error.txt');
        file_id = fopen(error_file, 'w');
        fprintf(file_id, '%s\n', getReport(this_error, 'extended', 'hyperlinks', 'off'));
        fclose(file_id);
    end
end