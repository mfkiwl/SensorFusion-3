function rosbag_read(folder_path)
% rosbag_read(folder_path) moves to the given folder path and reads all the
% rosbag files and perform data processing. Retrieved data are concatenated
% for convenience.
%
% Implemented by JinHwan Jeon, 2023
    
    % Access folder
    cd(folder_path); files = dir;
    n = length(files);

    % Read and process .bag files
    outputs = {};
    for i=1:n
        file_name = files(i).name;
        s = dir(file_name);
        size_ = s.bytes;
        if length(file_name) > 3
            if strcmp(file_name(end-2:end),'bag')
                disp(['Processing ',file_name, ': ',num2str(1e-6 * size_),' MB'])
                output = rosbag_convert(strcat(folder_path,file_name));
                outputs = [outputs {output}];
            end
        end
    end
    
    % Concatenate output from rosbag_convert
    output = catstruct(outputs);
    save('bag_data.mat','output');
end

