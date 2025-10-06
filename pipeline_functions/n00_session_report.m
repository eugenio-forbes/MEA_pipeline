%%% Function to generate a report of the timing of the
%%% pipeline processing of a recording session
s
function n00_session_report(root_directory, subject, folder, processing_times)

info_directory = fullfile(root_directory, 'MEA_database', subject, folder, 'info');
report_file = fullfile(info_directory, 'processing_times.txt');

if isfile(report_file)
    delete(report_file);
end

file_id = fopen(report_file, 'w')

fprintf(file_id, 'Processing times:\n');
fprintf(file_id, 'Info gathering and binary file creation:\n');
fprintf(file_id, 'Step 1: %.2f minutes\n', processing_times.step1);
fprintf(file_id, 'Combinato spike detection and clustering :\n');
fprintf(file_id, 'Step 2: %.2f minutes\n', processing_times.step2);
fprintf(file_id, 'Step 3: %.2f minutes\n', processing_times.step3);
fprintf(file_id, 'Step 4: %.2f minutes\n', processing_times.step4);
fprintf(file_id, 'Maxwell spike processing and clustering with combinato:\n');
fprintf(file_id, 'Step 2: %.2f minutes\n', processing_times.step2_maxwell_data);
fprintf(file_id, 'Step 3: %.2f minutes\n', processing_times.step3_maxwell_data);
fprintf(file_id, 'Step 4: %.2f minutes\n', processing_times.step4_maxwell_data);
fprintf(file_id, 'Kilosort processing and clustering:\n');
fprintf(file_id, 'Step 5: %.2f minutes\n', processing_times.step5);
fprintf(file_id, 'Processing of kilosort results:\n');
fprintf(file_id, 'Step 6: %.2f minutes\n', processing_times.step6);
fprintf(file_id, 'Cleanup:\n');
fprintf(file_id, 'Step 7: %.2f minutes\n', processing_times.step7);
fprintf(file_id, 'Total time: %.2f hours\n', processing_times.overall_time);

fclose(file_id);

end