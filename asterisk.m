function b = asterisk()
%adding asterisk to existing .dat files
data_dir_str= ['/Users/polinavanyukov/Box Sync/Project Trust Game/regs/13-Jun-2016/'];
data_dump_str = '/Users/polinavanyukov/Box Sync/Project Trust Game/regs/asterisked regs/';

if ~exist(data_dump_str,'file')
    mkdir(data_dump_str)
    fprintf('Creating specific reg folder in: %s\n\n',data_dump_str);
end

cd(data_dir_str)
files = dir('*.dat');
num_of_subjects = length(files);

parfor index = 1:num_of_subjects
    filename=files(index).name;
    fprintf('File processing: %s\n', filename);
    x = load(filename);
    block1=num2cell(x(1:48,:));
    block2=num2cell(x(49:96,:));
    block3=num2cell(x(97:144,:));
    block4=num2cell(x(145:192,:));
    ast = {'*', '*', '*'};
    c = [block1; ast; block2; ast; block3; ast; block4];
    %writetable(cell2table(c), [data_dump_str filename],'Delimiter','\t');
    dlmcell([data_dump_str filename],c,'\t');
end


return
