function make_review_combinato_files(varargin)
if isempty(varargin)
    rootDir = '/project/TIBIR/Lega_lab/shared/lega_ansir/Single_Units_BR';
else
    rootDir = varargin{1};
end
%List directories
manualDir = fullfile(rootDir,'micros_pipeline/process_files/n07_manual');

%Check if folder with combinato files to be able to open GUI exists
combinatoDir = fullfile(manualDir,'combinato_files');
if ~isfolder(combinatoDir)
    mkdir(combinatoDir);
end

%Create a file for fake microelectrode recording in the format used by combinato
fake_hdf5 = fullfile(combinatoDir,'NS6_001_fake.hdf5');
if ~isfile(fake_hdf5)
    fake_data = int16(sin(linspace(1,1800000,1800000))); %One minute
    fake_sr = 30000; %At sampling rate of 30kHz
    h5create(fake_hdf5,'/data',length(fake_data),'Datatype',class(fake_data));
    h5write(fake_hdf5,'/data',fake_data);
    h5create(fake_hdf5,'/sr',length(fake_sr),'Datatype',class(fake_sr));
    h5write(fake_hdf5,'/sr',fake_sr);
end

%Create folder that holds spike data
fake_data_dir = strrep(fake_hdf5,'.hdf5','');
if ~isfolder(fake_data_dir)
    mkdir(fake_data_dir);
end

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
%%%     -types                :        2 x n_groups   : (int16) for each group (1st row), the type (2nd row) (-1 artifact, 0 unclustered, 1 multi-unit, 2 single unit)
%%%     -types_orig           :        2 x n_groups   : (int16) same as above, before manual clustering (not used)
%%%     -artifacts            :        2 x n_classes  : (int64) for each class (1st row), whether it is artifact (1) or not (0) (2nd row)
%%%     -artifacts_prematch   :        2 x n_classes  : (int64) same as above, before manual clustering (not used)
%%%     -groups               :        2 x n_classes  : (int16) for each class (1st row), which group (2nd row)
%%%     -groups_orig          :        2 x n_classes  : (int16) same as above, before manual clustering (not used)
%%%     -classes              : (n_spikes-skipped) x 1: (uint16) class number for each spike
%%%     -distance             : (n_spikes-skipped) x 1: (single) template matching distances
%%%     -index                : (n_spikes-skipped) x 1: (uint32) python index (0:length(times)-1) of the spike times corresponding to spike (some skipped)
%%%     -matches              : (n_spikes-skipped) x 1: (int8) 0 (SPC), 1 (1st template matching), 2 (second template matching)


%Create file with fake spike data
fake_spike_file = fullfile(fake_data_dir,'data_NS6_001_fake.h5');
if ~isfile(fake_spike_file)
    fake_spikes = single(repmat(sin(linspace(1,60,60))',1,150));
    fake_times = double(linspace(1,60000,150))';
    fake_artifacts = int8(zeros(150,1));
    
    h5create(fake_spike_file,'/neg/spikes',size(fake_spikes),'Datatype',class(fake_spikes));
    h5write(fake_spike_file,'/neg/spikes',fake_spikes);
    h5create(fake_spike_file,'/neg/times',length(fake_times),'Datatype',class(fake_times));
    h5write(fake_spike_file,'/neg/times',fake_times);
    h5create(fake_spike_file,'/neg/artifacts',length(fake_artifacts),'Datatype',class(fake_artifacts));
    h5write(fake_spike_file,'/neg/artifacts',fake_artifacts);
    
    h5create(fake_spike_file,'/pos/spikes',size(fake_spikes),'Datatype',class(fake_spikes));
    h5write(fake_spike_file,'/pos/spikes',fake_spikes);
    h5create(fake_spike_file,'/pos/times',length(fake_times),'Datatype',class(fake_times));
    h5write(fake_spike_file,'/pos/times',fake_times);
    h5create(fake_spike_file,'/pos/artifacts',length(fake_artifacts),'Datatype',class(fake_artifacts));
    h5write(fake_spike_file,'/pos/artifacts',fake_artifacts);
end

%Create folder with fake sorting results for those spikes
fake_sort_dir = fullfile(fake_data_dir,'sort_neg_fake');
if ~isfolder(fake_sort_dir)
    mkdir(fake_sort_dir);
end

%Create fake sorting results file
fake_sort_file = fullfile(fake_data_dir,'sort_cat.h5');
if ~isfile(fake_sort_file)
    fake_types = int16([0,1;-1,1]);
    fake_artifacts = int64([0,1;1,0]);
    fake_groups = int16([0,1;0,1]);
    fake_classes = uint16(ones(150,1));
    fake_distance = single(ones(150,1));
    fake_index = uint32(1:150);
    fake_matches = int8(zeros(150,1));
    
    h5create(fake_sort_file,'/types',size(fake_types),'Datatype',class(fake_types));
    h5write(fake_sort_file,'/types',fake_types);
    h5create(fake_sort_file,'/types_orig',size(fake_types),'Datatype',class(fake_types));
    h5write(fake_sort_file,'/types_orig',fake_types);
    h5create(fake_sort_file,'/artifacts',size(fake_artifacts),'Datatype',class(fake_artifacts));
    h5write(fake_sort_file,'/artifacts',fake_artifacts);
    h5create(fake_sort_file,'/artifacts_prematch',size(fake_artifacts),'Datatype',class(fake_artifacts));
    h5write(fake_sort_file,'/artifacts_prematch',fake_artifacts);
    h5create(fake_sort_file,'/groups',size(fake_groups),'Datatype',class(fake_groups));
    h5write(fake_sort_file,'/groups',fake_groups);
    h5create(fake_sort_file,'/groups_orig',size(fake_groups),'Datatype',class(fake_groups));
    h5write(fake_sort_file,'/groups_orig',fake_groups);
    h5create(fake_sort_file,'/classes',length(fake_classes),'Datatype',class(fake_classes));
    h5write(fake_sort_file,'/classes',fake_classes);
    h5create(fake_sort_file,'/distance',length(fake_distance),'Datatype',class(fake_distance));
    h5write(fake_sort_file,'/distance',fake_distance);
    h5create(fake_sort_file,'/index',length(fake_index),'Datatype',class(fake_index));
    h5write(fake_sort_file,'/index',fake_index);
    h5create(fake_sort_file,'/matches',length(fake_matches),'Datatype',class(fake_matches));
    h5write(fake_sort_file,'/matches',fake_matches);
end
end
    