

%main> 
function muneshwardatanalysisfinal
%file_name1,loc_data_lfp,TD_col,TD_row,sum_col,sum_rows,total_RL,total_LR,sum_DG_RL,sum_DG_LR,lfp_res
current_dir=pwd;
curr_dir='E:\_ANN\data\ephysdata';
dirinfo = dir(curr_dir);
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
tf = ismember( {dirinfo.name}, {'.', '..'});
dirinfo(tf) = [];

for folder_no=1 %20:20 %length(dirinfo)
    datadate=dirinfo(folder_no).name;
    part_path_name = 'E:\_ANN\data\ephysdata';
    data_path=strcat(part_path_name,datadate,'\');
    list = dir(fullfile(data_path, '*.mat')); % f = fullfile(filepart1,...,filepartN) builds a full file specification from the specified folder and file names.
    file_name1 = {list.name};
    str  = sprintf('%s#', file_name1{:});
    num=sscanf(str, '%*5c_%*4c_%*u-%*3c-%*u_%*u_%*u_%*u-%u.mat#');
    [dummy, index] = sort(num);
    file_name1 = file_name1(index);
    % curr_dir=pwd;
    
    
    
    cd(data_path);
    cd(current_dir);
    %eval(sprintf('cd %s',data_path));
    %eval(sprintf('cd %s',current_dir));
    
    for k =1:length(file_name1)
        file_name = file_name1{k};
        fprintf('\nNow analyzing %s\n', file_name);
        muneshwar(data_path,file_name);
        %         lfp_analysis(data_path,file_name)
        % class_unit(data_path,file_name)
        eval(sprintf('cd %s',current_dir));
    end
end
end