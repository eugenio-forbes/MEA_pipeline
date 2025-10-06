%%% MEA Pipeline, Step 3: Use Combinato to Cluster Spike Waveforms
%%%
%%% This function will use presaved spike clustering parameters
%%% and use system commands to execute python-based Combinato
%%% to cluster spike waveforms for every channel.

function [error_flag, error_message] = n03_combinato_cluster(varargin)
if isempty(varargin)
    root_directory = '/path/to/MEA_pipeline/parent_directory';
    subject = 'SC000';
    folder = 'yyyy-mm-dd_network';
    is_maxwell_data = true; 
else
    root_directory = varargin{1};  %%% (character array) Parent directory containing /MEA_pipeline and /MEA_database.
    subject = varargin{2};         %%% (character array) Subject code in 'SC000' format.
    folder = varargin{3};          %%% (character array) Folder name with data in 'yyyy-mm-dd_scan-name#' format. Raw recording in /folder/raw/data.raw.h5
    is_maxwell_data = varargin{4}; %%% (true/false)      Whether spike waveforms were extracted based on MW spike times.
end

%%% List directories and files
data_directory = fullfile(root_directory, 'MEA_database', subject, folder);
combinato_directory = fullfile(data_directory, 'combinato');
template_directory = fullfile(root_directory, 'MEA_pipeline/base_code/combinato_templates');

if is_maxwell_data
    combinato_directory = [combinato_directory, '_maxwell_data'];
    cluster_template_file = fullfile(template_directory, 'template_cluster_maxwell_data.txt');
else
    cluster_template_file = fullfile(template_directory, 'template_cluster.txt');
end

cluster_template = fileread(cluster_template_file);
options_text_file =  fullfile(template_directory, 'custom_options.txt');
options_text = fileread(options_text_file);

delete_previous_attempts(combinato_directory);    %%% Because duplicate jobs cause errors

%%% Print custom options to local_options.py file and make executable
options_file = fullfile(combinato_directory, 'local_options.py');
file_id = fopen(options_file, 'w'); 
fprintf(file_id, options_text); 
fclose(file_id);
system(sprintf('chmod +x %s', options_file));

%%% Values to be returned by function
error_flag = 0;
error_message = '';

%%% Entries for clustering template
entries = [];
entries = [entries; {combinato_directory}]; %%% combinato_directory first
label = 'options00';
entries = [entries; repelem({label}, 6, 1)]; %%% then 6 times label

%%% Print bash file for clustering spikes and make executable
cluster_text = sprintf(cluster_template, entries{:}); %%% Input entries to template
cluster_file = fullfile(combinato_directory, 'cluster.sh');
file_id = fopen(cluster_file, 'w');
fprintf(file_id, cluster_text);
fclose(file_id);
system(sprintf('chmod +x %s', cluster_file));

%%% Run clustering in system shell
[~, command_output] = system(sprintf('bash %s', cluster_file));

%%% End function execution if combinato outputs traceback or error
%%% Return combinato output and flag error
if contains(command_output, {'traceback', 'error'}, 'IgnoreCase', true)
    error_flag = 1;
    error_message = command_output;
end

delete(cluster_file);
delete(options_file);
rmdir(fullfile(combinato_directory, '__pycache__'), 's');

end


function delete_previous_attempts(combinato_directory)
    %%%Delete previous attempt files if there are any. Repeat run would mess things up.
    
    %%% First delete overview and pycache folders
    pycache_directory = fullfile(combinato_directory, '__pycache__');
    if isfolder(pycache_directory)
        rmdir(pycache_directory, 's');
    end
    
    overview_directory = fullfile(combinato_directory, 'overview');
    if isfolder(overview_directory)
        rmdir(overview_directory, 's');
    end
    
    %%% List remaining folders and files
    invalid = {'.', '..', 'do_sort_neg.txt', 'do_sort_pos.txt'};
    directory_list = dir(combinato_directory);
    directory_list = directory_list(~ismember({directory_list.name}, invalid));
    
    %%% Anything that isn't a folder in the first level can be deleted
    %%% Except do_sort_neg.txt and do_sort_pos.txt files
    combinato_files = directory_list(~[directory_list.isdir]);
    combinato_files = fullfile({combinato_files.folder}, {combinato_files.name});
    for idx = 1:length(combinato_files)
        this_file = combinato_files{idx};
        delete(this_file);
    end
    
    %%% Loop through combinato directories for each file and delete log.txt
    %%% and all directories within it, thus only keeping file with spikes.
    directories = directory_list([directory_list.isdir]);
    combinato_folders = fullfile({directories.folder}, {directories.name});
    for idx = 1:length(combinato_folders)
        this_folder = combinato_folders{idx};
        log_file = fullfile(this_folder, 'log.txt');
        if isfile(log_file)
            delete(log_file);
        end
        subfolders = dir(this_folder);
        subfolders = subfolders(~ismember({subfolders.name}, invalid)); 
        subfolders = subfolders([subfolders.isdir]);
        subfolders = fullfile({subfolders.folder}, {subfolders.name});
        for jdx = 1:length(subfolders)
            rmdir(subfolders{jdx}, 's');
        end 
    end
end