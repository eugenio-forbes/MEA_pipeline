function update_MEAs(varargin)
if isempty(varargin)
    root_directory = '/path/to/MEA_pipeline/parent_directory';
    do_combinato = true;
    do_combinato_maxwell = true;
    do_kilosort = true;
    do_plots = false;
else
    root_directory = varargin{1};                 %%% Parent directory with pipeline and database folders (string/character array).
    do_combinato = logical(varargin{2});          %%% Whether to detect spikes and sort with combinato methods. (true/false)
    do_combinato_maxwell = logical(varargin{3});  %%% Whether to sort spikes detected by Maxwell software with combinato. (true/false)
    do_kilosort = logical(varargin{4});           %%% Whether to use kilosort3 for detection and sorting. (true/false)
    do_plots = logical(varargin{5});              %%% Whether to make plots at each step used. Slower processing. (true/false)
    %%% Can only do plots using matlab/2023b+ because of colormap selection.
end

%%% Declare directories
code_directory = fullfile(root_directory, 'MEA_pipeline');
addpath(genpath(code_directory));
data_directory = fullfile(root_directory, 'MEA_database');
error_directory = fullfile(data_directory, 'error_logs');

%%% List all Maxlab raw file paths
raw_file_paths = dir(fullfile(data_directory, '*/*/raw/*.raw.h5'));

if isempty(raw_file_paths)
    error('Did not find any raw files in %s.\n', data_directory);
end

raw_file_paths = {raw_file_paths.folder};
raw_file_paths = strrep(raw_file_paths, [data_directory '/'], '');
raw_file_paths = strrep(raw_file_paths, '/raw', '')';

%%% Check for existence of table listing all sessions and if it does exist
%%% exclude from paths those that have already been completely processed

session_list = fullfile(data_directory, 'all_sessions.mat');
if isfile(session_list)
    load('session_list', 'all_sessions');
    all_sessions = all_sessions(all_sessions.complete, :);
    completed_paths = fullfile(all_sessions.subject, all_sessions.folder);
    raw_file_paths = raw_file_paths(~ismember(raw_file_paths, completed_paths));
end

%%% Gather subject and folder names for all paths
subjects = cellfun(@(x) {x(1:strfind(x, '/') - 1)}, raw_file_paths);
folders = cellfun(@(x) {x(strfind(x, '/') + 1:end)}, raw_file_paths);

%%% Loop through all folders to process through pipeline
timer_overall = tic;
error_count = 0;
for idx = 1:length(folders)
    subject = subjects{idx};
    folder = folders{idx};
    try
    
        processing_times = struct;
        timer_session = tic;

        %%% Step 1. All steps will use output of this function: All maxwell
        %%% raw files have information gathered, channels are mapped and
        %%% given a row ID for input of individual channel results into a
        %%% column vector so that they can be reshaped into the shape of the
        %%% MEA and plotted with imagesc/contourf. The signal is demeaned
        %%% and saved as int16 binary file mainly for use by kilosort and by other
        %%% steps.
        
        timer_1 = tic;
        n01_maxwell2kilosort(root_directory, subject, folder, do_plots);
        processing_times.step1 = toc(timer_1)/60; clear timer_1
        
        %%% Steps 2-4. Two different options that can both be carried
        %%% out based on variables do_combinato and do_combinato_maxwell.   
        %%% 1) Detect and sort spikes with combinato.
        %%% 2) Extract spikes based on Maxwell software spike times
        %%%    and methods, and sort with combinato.
        %%% - Channel/neuron/class clustering quality metrics gathered. 
        %%% - Spikes saved in combinato format for access.
        
        if do_combinato
            is_maxwell_data = false;
            
            %%% Step 2. Spike Extraction. Positive and negative spikes.
            timer_2 = tic;
            n02_extract_spikes_combinato(root_directory, subject, folder);
            processing_times.step2 = toc(timer_2) / 60; clear timer_2
            
            %%% Step3. Combinato clustering
            timer_3 = tic;
            n03_combinato_cluster(root_directory, subject, folder, is_maxwell_data);
            processing_times.step3 = toc(timer_3) / 60; clear timer_3
            
            %%% Step4. Channel/neurons/classes metrics gathered.
            timer_4 = tic;
            n04_combinato_data(root_directory, subject, folder, is_maxwell_data);
            processing_times.step4 = toc(timer_4) / 60; clear timer_4
            
            if do_plots
                n00_plot_combinato_metrics(root_directory, subject, folder, is_maxwell_data);
            end
        else
            processing_times.step2 = 0;
            processing_times.step3 = 0;
            processing_times.step4 = 0;
        end
        
        if do_combinato_maxwell
            is_maxwell_data = true;
            
            %%% Step 2. Spike Extraction. Only negative spikes. Only highpass filtered over 300 Hz.
            timer_2 = tic;
            n02_extract_spikes_combinato_maxwell(root_directory, subject, folder);
            processing_times.step2_maxwell = toc(timer_2) / 60; clear timer_2
            
            %%% Step3. Combinato clustering
            timer_3 = tic;
            n03_combinato_cluster(root_directory, subject, folder, is_maxwell_data);
            processing_times.step3_maxwell = toc(timer_3) / 60; clear timer_3
            
            %%% Step4. Channel/neurons/classes metrics gathered.
            timer_4 = tic;
            n04_combinato_data(root_directory, subject, folder, is_maxwell_data);
            processing_times.step4_maxwell = toc(timer_4) / 60; clear timer_4
            
            %%% Plot channel/neurons/classes metrics
            if do_plots
                n00_plot_combinato_metrics(root_directory, subject, folder, is_maxwell_data)
            end
        else
            processing_times.step2_maxwell = 0;
            processing_times.step3_maxwell = 0;
            processing_times.step4_maxwell = 0;
        end
        
        %%% Step 5 and 6. Process raw binary file to detect and cluster
        %%% spikes with kilosort. Gather results to create files for each
        %%% channel that contained units and spike time series for each unit.
        if do_kilosort
            timer_5 = tic;
            n05_kilosort_cluster(root_directory, subject, folder);
            processing_times.step5 = toc(timer_5)/60; clear timer_5
            timer_6 = tic;
%             n06_get_kilosort_units(root_directory, subject, folder);
            processing_times.step6 = toc(timer_6)/60; clear timer_6
        else
            processing_times.step5 = 0;
            processing_times.step6 = 0;
        end
        
        %%% Step 7. Delete all intermediate files (raw binary)
        timer_7 = tic;
        n07_clean_up(root_directory, subject, folder);
        processing_times.step7 = toc(timer_7); clear timer_7
        
        processing_times.overall_time = toc(timer_session)/3600; clear timer_session
        
        %%% Make report for processing of this session
        n00_session_report(root_directory, subject, folder, processing_times);
        
    catch this_error
        error_count = error_count + 1;
        this_error_directory = fullfile(error_directory, subject, folder);
        if ~isfolder(this_error_directory)
            mkdir(this_error_directory);
        end
        error_file = fullfile(this_error_directory, 'pipeline_error.txt');
        file_id = fopen(error_file, 'w');
        fprintf(file_id, '%s\n', getReport(this_error, 'extended', 'hyperlinks', 'off'));
        fclose(file_id);
    end    
end

%%% Gather all session_info, recording_list files, and make a single file for next time.
n08_get_all_session_info(root_directory);

running_time = toc(timer_overall);
fprintf('Processed %d incomplete sessions in %.02f hours.\n', length(folders), running_time);
fprintf('Out of the %d sessions processed, %d had an error.\n', length(folders), error_count);

end