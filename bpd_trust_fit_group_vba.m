function bpd_trust_fit_group_vba
jj=1;
hh=1;

%Load the sceptic config files and initialize the tracking data.
%As a new user you will have to create the config files (see
%https://github.com/DecisionNeurosciencePsychopathology/temporal_instrumental_agent
%for more help) & set the paths to said file.
task_data=initialize_task_tracking_data('trust_bpd');


%Create hallquist overrride?

%Quick username check, and path setting, this may have to change depending
%on the machine you are currently working on!
os = computer;
if strcmp(os(1:end-2),'PCWIN')
    %     %data_dir = 'E:/Users/wilsonj3/Google Drive/skinner/data - IP/eprime/bpd_trust/';
    data_dir='subjects';
    %files=glob('subjects/*/*_scan_*.txt');
    %     %files = glob([data_dir '/*/*scan_*.txt']);
    local_data_dir = 'E:\data\bpd_trust\vba_out'; %Where the vba workspace will be saved
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

%Just for now! 10/10/2016  Takes to long already processed!
% for i = 1:length(files)
%     [~,ids(i)]=bpd_trustbehavior(data_dir,files{i});
% end

%save trust_bpd_ids ids

%JON! remove the loading of the ids, grab the ids from the subjects dir!!!


%destination folder on Thorndike (server) to house regressors
dest_folder='/Volumes/bek/bsocial/bpd_trust/regs'; 

%get file paths
scriptName = mfilename('fullpath');
[currentpath, ~, ~]= fileparts(scriptName);

%Get subject ids
ids=dir('subjects');
ids=regexp({ids.name},'\d+','match');
ids=cellfun(@str2double,[ids{:}]);

%load('trust_bpd_ids.mat') %This needs to be automatically generated and or updated!
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
model = 'subjectCounterfactual';

%% main loop
L = [];
for i = 1:length(ids)
    %     filename=files{i};
    %     fprintf('File processing: %s\n', filename);
    id = ids(i);
    datalocation = [data_dir '/' num2str(id)];
    
    %Update task_tracking data
    task_data.behave_completed=1;
    
    %only process behavior for those not processed yet
    if ~exist([datalocation '/' sprintf('trust%d.mat',id)])
        file = glob([datalocation  '/' '*_scan_*.txt']);
        bpd_trustbehavior(data_dir,file{:});
    end
    
    %Update task_tracking data 
    task_data.behave_processed=1;
    
    %datalocation = [data_dir];
    %[posterior, out] = bpd_trust_Qlearning_ushifted(id, counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret);
    [~, out] = bpd_trust_Qlearning_ushifted(id, datalocation,...
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
        moveregs(currentpath,num2str(id),dest_folder);
        
        %write the task data to file
        record_subj_to_file(id,task_data)
        
    catch exception
        %write the task data to file
        record_subj_to_file(id,task_data)
        
        %Error catcher
        errorlog('bpdtrust',id,exception)
        
    end
end
    
    %Need to add modified PE scripts, potentially need to add model loop
    L_name = sprintf('test_L_model%s_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d_regret%d',model, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret);
    save(char(L_name), 'L'); %Just saveing L's
