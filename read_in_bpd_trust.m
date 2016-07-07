function b = read_in_bpd_trust(fname, procname, fields, startoffset, endoffset)
 %fields = {'ConditionOrder','identity','Running','TrusteeDecides','TrialNumber','OrderRS','exchangeNum','Reversal', 'partnerchoice.OnsetTime','partnerchoice.RTTime','partnerchoice.RT','PartDecides'};
 %startoffset = 14; %from TrialProc occurrence to the beginning of the block containing relevant info
 %endoffset = 28;    %from TrialProc occurrence to the end of the block containing relevant info
if (~exist(fname,'file')) 
    error('MATLAB:eprimeread:FileNotFound','File:\n\t%s\n was not found.\n',fname); 
else
	b.fname = fname; 
end
% check for file format
if(checkFileBOM(fname)) % file format is UTF-16 IEEE-LE
    fprintf('\tfile: ''%s''\n\tis in UTF16-LE format\n',fname);
    fprintf('\t\tconverting to ASCII...');
    utf2ascii(b.fname,'fast'); % convert it (fast = no need to verify BOM) %was uni2ascii()
    fprintf('\t done\n');
end
%read the raw text file
fid = fopen(b.fname,'r');

if(fid < 0)
    warning('eprimeread:open_file','Unable to open file ''%s'' ...moving on',fname);
    fclose(fid);
    return
else
    rawtxt = textscan(fid,'%s','delimiter',':\n','whitespace',''); % why '' for whitespace?
    rawtxt = rawtxt{:};
    fclose(fid);
end
% strip out non-fields
q_keep = ( cellfun(@numel,rawtxt) > 0 ); % set empty char. entries to 'false'
q_keep = ( q_keep & cellfun(@isempty,strfind(rawtxt,'***')) );
txt = cellfun(@strtrim,rawtxt(q_keep),'UniformOutput',false); % txt is a cell array where each row is a single string from the raw data file

% get lines on which the correct procedure occurs
trials=0; % trials variable final value is the number of times the indicated procedure (TrialProc) occurs in the raw data file
proclines=0; % proclines is an array of lines in the raw data file on which the indicated procedure (TrialProc) occurs

for ct=1:length(txt)-1 
    if (strcmp(procname,txt{ct}))         
        trials = trials + 1; 
		proclines(trials) = ct;
    end
end

data = cell(trials, length(fields));
numdata = zeros(trials, length(fields))-999;
line = 0;

% Reading in all of the data
% note that if a partial match of the field name is needed, use strfind
% instead of strcmp below
for ind = 1:trials
    startline = proclines(ind)-startoffset*2;
    endline = proclines(ind)+endoffset*2;
    for line=startline:endline
        for field=1:length(fields)
        if(strcmp(strtrim(char(txt(line))),char(fields(field))))
            data(ind,field)=txt(line+1);
            numval=str2double(char(data(ind,field)));
            if isnan(numval) 
                numval=[]; 
            elseif (~isempty(numval))
                numdata(ind,field)=numval;			
            else
	  			numdata(ind,field)=-999;
            end
        end
        end             
    end
    
end



if(isempty(data{end,1}))
	data=data(1:end-1,:);
	numdata=numdata(1:end-1,:);
end

for ct=1:length(fields)
    if max(numdata(:,ct))==-999
    	curdata=data(:,ct);
  	else
    	curdata=numdata(:,ct);
  	end
  	fldname=strrep(char(fields(ct)),'.','_');
   %fldname=strrep(char(fields(ct)),':','');
   % fldname=genvarname(char(fields(ct)));
  	b=setfield(b,fldname,curdata);
end

    
return

    
            