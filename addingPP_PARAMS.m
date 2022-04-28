
function addingPP_PARAMS(date)
cd('C:\Users\àdmin\Desktop\adarsh\MATLAB\sorted\rss_ssa\')

Experiment_date=date;
datadate=Experiment_date;
current_dir=pwd;

new_path=strcat('C:\Users\àdmin\Desktop\adarsh\MATLAB\sorted\rss_ssa\',datadate);
o_path=strcat('C:\Users\àdmin\Desktop\adarsh\MATLAB\data\rss_ssa\',datadate);
cd(o_path)
files=dir('*.mat');
for kk=1:length(files)
    file_name = files(kk).name;
    vars = whos('-file',file_name);
    load(file_name,vars(1).name)
    cd(new_path);
    add_unit=strcat(file_name(1:end-4),'_unit_record.mat');
    load(add_unit);
    save(add_unit,'PP_PARAMS','unit_record','unit_record_spike')%,'spikeyn')
    clear PP_PARAMS unit_record unit_record_spike
    cd(new_path)
    cd(o_path)
    kk
end

end