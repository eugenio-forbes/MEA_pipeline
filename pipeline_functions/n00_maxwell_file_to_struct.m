%%% The following three functions loop recursively through the hierarchical
%%% data file groups, to load datasets into a hierarchical matlab structure

function maxwell_struct = maxwell_file_to_struct(maxwell_file)

maxwell_struct = struct;

maxwell_file_info = h5info(maxwell_file);

start_path = '/';
groups = maxwell_file_info.Groups;
datasets = maxwell_file_info.Datasets;

maxwell_struct = add_datasets(maxwell_file, maxwell_struct, start_path, datasets);

maxwell_struct = add_groups(maxwell_file, maxwell_struct, start_path, groups);

end


function parent_struct = add_groups(maxwell_file, parent_struct, parent_path, groups)
for idx = 1:length(groups)

    group_name = groups(idx).Name;
    delimiter = strfind(group_name, '/');
    
    if ~isempty(delimiter)
        group_name = group_name(delimiter(end) + 1:end);
    end
    
    if ~strcmp(group_name, 'spikes')
        group_path = fullfile(parent_path, group_name);
        
        group_datasets = groups(idx).Datasets;
        group_groups = groups(idx).Groups;
        
        group_struct = struct;
        group_struct = add_datasets(maxwell_file, group_struct, group_path, group_datasets);
        group_struct = add_groups(maxwell_file, group_struct, group_path, group_groups);
        
        parent_struct.(group_name) = group_struct;
    end
    
end
end


function parent_struct = add_datasets(maxwell_file, parent_struct, parent_path, datasets)
for idx = 1:length(datasets)

    dataset_name = datasets(idx).Name;
    dataset_path = fullfile(parent_path, dataset_name);
    dataset_name = strrep(dataset_name, ' ', '_');
    
    try
        if strcmp(dataset_name, 'raw')
            dataset = dataset_path;
        else
            dataset = h5read(maxwell_file, dataset_path);
            if iscell(dataset) && numel(dataset) == 1
                if isnumeric(dataset{1})
                    dataset = dataset{1};
                end
            end
        end
    catch
        dataset = [];
    end
    
    parent_struct.(dataset_name) = dataset;
    
end
end