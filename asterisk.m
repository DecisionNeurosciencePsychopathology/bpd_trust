function asterisk(id)
%adding asterisk to existing .dat files
data_dir_str= 'regs/';
data_dump_str = sprintf('asterisked regs/%d/',id);

cd(data_dir_str)

if ~exist(data_dump_str,'file')
    mkdir(data_dump_str)
    fprintf('Creating id specific asterisked reg folder in: %s\n\n',data_dump_str);
end

%/Users/polinavanyukov/Box Sync/Project Trust Game/
% if ~exist(data_dump_str,'file')
%     mkdir(data_dump_str)
%     fprintf('Creating specific reg folder in: %s\n\n',data_dump_str);
% end

cd(sprintf('%d', id))
files = dir(sprintf('%d*.dat', id));
num_of_subjects = length(files);


for index = 1:num_of_subjects
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
    cd ..
    data_dump_str = sprintf('asterisked regs/%d/',id);
    dlmcell([data_dump_str '/' filename],c,'\t');
    cd(sprintf('%d', id))
end

% Sometimes it need to be cd ../..
cd ..


return

