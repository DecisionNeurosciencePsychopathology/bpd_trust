ids_with_four_runs = regexp(glob('C:\kod\bpd_trust\regs\asterisked regs\*'),'\d+','match');

for i = 1:length(ids_with_four_runs)
   id=ids_with_four_runs{i};
   asterisk(str2double(id{:})); 
   %Copy censor reg
   copyfile(['C:\kod\bpd_trust\regs\' id{:} filesep id{:} 'to_censor_motion_corr'],['C:\kod\bpd_trust\regs\asterisked regs\' id{:} filesep id{:} 'to_censor_motion_corr']);
end