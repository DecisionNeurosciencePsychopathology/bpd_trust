function bpd_trust_fit_group_vba
 jj=1;
 hh=1;

 


%Create hallquist overrride?

%Quick username check, and path setting, this may have to change depending
%on the machine you are currently working on!
os = computer;
if strcmp(os(1:end-2),'PCWIN')
%     %data_dir = 'E:/Users/wilsonj3/Google Drive/skinner/data - IP/eprime/bpd_trust/';
     data_dir='subjects';
      %files=glob('subjects/*/*_scan_*.txt');
%     %files = glob([data_dir '/*/*scan_*.txt']);
     local_data_dir = 'E:\data\bdp_trust\vba_out'; %Where the vba workspace will be saved
%     %local_data_dir = ''; %Use this options to not save the vba data
%     %processed_files = glob('E:/Users/wilsonj3/Google Drive/skinner/data - IP/eprime/bpd_trust/*/*.mat');
    
    
    
%     %data_dir = 'E:/Users/wilsonj3/Google Drive/skinner/data - IP/eprime/bpd_trust/';
%     data_dir='C:\Users\wilsonj3\Desktop\hallquist_trust';
% %     files=glob('subjects/*/*_scan_*.txt');
%     %files = glob([data_dir '/*/*scan_*.txt']);
%     %local_data_dir = 'E:\data\trust'
%     local_data_dir = 'E:\data\hallquist trust'; %Where the vba workspace will be saved -- hallquist
%     %local_data_dir = ''; %Use this options to not save the vba data
%     %processed_files = glob('E:/Users/wilsonj3/Google Drive/skinner/data - IP/eprime/bpd_trust/*/*.mat');
    
else
    [~, me] = system('whoami');
    me = strtrim(me);
    if strcmp(me,'polinavanyukov')==1
        data_dir = '';
        files = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/beha_behavior/');
        local_data_dir = '';
    else
        files = glob('?');
      

    end
end

%Note is reprocessing the files over and over again takes to long implment
%the idea that we'll only process the new ones
%if length(raw_files) ~= length(processed_files)
%    for 


%Just for now! 10/10/2016  Takes to long already processed!
% for i = 1:length(files)
%     [~,ids(i)]=bpd_trustbehavior(data_dir,files{i});
% end

%save trust_bpd_ids ids
load('trust_bpd_ids.mat') %This needs to be automatically generated and or updated!
%load('hallquist_trust_ids.mat')
% % % % %% chose models to fit
% % % % modelnames = {'ushifted_bpd_trust_Qlearning'};

%% set parameters
% nbasis = 4;
% multinomial = 1;
%counter = 1;                    %using counterfactual feedback
multisession = 0;               %modelling runs separately
fixed_params_across_runs = 1;   
sigma_kappa = 1;                %kappa (or action bias) parameter
reputation_sensitive = 0;       %modelling trustees' reputation
humanity = 0;                   %modelling humanity
valence_p = 0;                  %modelling valence of trustees
valence_n = 0;                  
assymetry_choices = 0;          %modelling assymetry in choices
regret = 0;
model = 'null';

%% main loop
L = [];
for i = 1:length(ids)
%     filename=files{i};
%     fprintf('File processing: %s\n', filename);
    id = ids(i);
    datalocation = [data_dir '/' num2str(id)];
    
    if ~exist([datalocation '/' sprintf('trust%d.mat',id)])
        file = glob([datalocation  '/' '*_scan_*.txt']);
        bpd_trustbehavior(data_dir,file{:});
    end
    
    %datalocation = [data_dir];
    %[posterior, out] = bpd_trust_Qlearning_ushifted(id, counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret);
    [posterior, out] = bpd_trust_Qlearning_ushifted(id, datalocation,...
        'multisession', multisession, 'fixed', fixed_params_across_runs,...
        'sigmakappa', sigma_kappa, 'reputation_sensitive', reputation_sensitive,...
        'humanity', humanity, 'valence_p', valence_p, 'valence_n', valence_n,...
        'assymetry_choices', assymetry_choices, 'regret', regret, ...
        'local_data_dir', local_data_dir, 'model', model);
    L(i) = out.F;
    
try
    %NB: IF we introduce a model loop this should only be ran once, though the PE regs will have to be ran more than once
    b=bpd_trustmakeregressor_group(id); 
    
    %move the regressor files to thorndike
    newfolder='/Volumes/bek/bsocial/bpd_trust/regs'; %folder to be place in within thorndike
    
    %get file paths
    scriptName = mfilename('fullpath');
    [currentpath, filename, fileextension]= fileparts(scriptName);
    moveregs(currentpath,num2str(id),newfolder);
    
    %write the ids that successfully ran into a cell
    ID(jj,1)=id;
    
    
    task={'bpd_trust'};
    Task{jj,1}=task; 
    
    trialdone=fopen('idlog_bpdtrust.txt', 'a+');
    trialdone=fscanf(trialdone,'%d');
    
    trialdone1=0;
    
    for aa=1:length(trialdone)
        if trialdone(aa,1) == id
            trialdone1=1;
        end
    end
    
      
    if trialdone1 == 1
        td={'yes'};
    else
        td={'no'};
    end
    fMRI_Preprocess_Complete{jj,1}=td; 
    
    jj=jj+1;
    
    %turn completed cell into table
    tt=table(ID,Task,fMRI_Preprocess_Complete);
    save('completed','tt');
    
catch exception
    
    errorlog('bpdtrust',id,exception)
        
    %put IDs that didn't run into table
        ID2(hh,1)=id; 
    
        task={'bpd_trust'};
        Task2{hh,1}=task; 
        
        hh=hh+1;
        
        tt2=table(ID2,Task2);
        save('unable_to_run','tt2')
end


if exist('tt2')==0
    ID2=0;
    Task2={'bpd_trust'};
    tt2=table(ID2,Task2);
    save('unable_to_run','tt2')
end
    
end

%Need to add modified PE scripts, potentially need to add model loop



L_name = sprintf('test_L_model%s_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d_regret%d',model, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret);
save(char(L_name), 'L'); %Just saveing L's
