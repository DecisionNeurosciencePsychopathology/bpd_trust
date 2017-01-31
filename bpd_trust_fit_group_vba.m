function bpd_trust_fit_group_vba


%Quick username check, and path setting, this may have to change depending
%on the machine you are currently working on!
os = computer;
if strcmp(os(1:end-2),'PCWIN')
%     %data_dir = 'E:/Users/wilsonj3/Google Drive/skinner/data - IP/eprime/bpd_trust/';
%     data_dir='subjects';
%     files=glob('subjects/*/*_scan_*.txt');
%     %files = glob([data_dir '/*/*scan_*.txt']);
%     local_data_dir = 'E:\data\bdp_trust\vba_out'; %Where the vba workspace will be saved
%     %local_data_dir = ''; %Use this options to not save the vba data
%     %processed_files = glob('E:/Users/wilsonj3/Google Drive/skinner/data - IP/eprime/bpd_trust/*/*.mat');
    
    
    
    %data_dir = 'E:/Users/wilsonj3/Google Drive/skinner/data - IP/eprime/bpd_trust/';
    data_dir='C:\Users\wilsonj3\Desktop\hallquist_trust';
%     files=glob('subjects/*/*_scan_*.txt');
    %files = glob([data_dir '/*/*scan_*.txt']);
    %local_data_dir = 'E:\data\trust'
    local_data_dir = 'E:\data\hallquist trust'; %Where the vba workspace will be saved -- hallquist
    %local_data_dir = ''; %Use this options to not save the vba data
    %processed_files = glob('E:/Users/wilsonj3/Google Drive/skinner/data - IP/eprime/bpd_trust/*/*.mat');
    
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
%%%load('trust_bpd_ids.mat')
load('hallquist_trust_ids.mat')
%% chose models to fit
modelnames = {'ushifted_bpd_trust_Qlearning'};

%% set parameters
% nbasis = 4;
% multinomial = 1;
counter = 1;                    %using counterfactual feedback
multisession = 0;               %modelling runs separately
fixed_params_across_runs = 1;   
sigma_kappa = 1;                %kappa (or action bias) parameter
reputation_sensitive = 0;       %modelling trustees' reputation
humanity = 0;                   %modelling humanity
valence_p = 0;                  %modelling valence of trustees
valence_n = 0;                  
assymetry_choices = 0;          %modelling assymetry in choices
regret = 0;

%% main loop
L = [];
for i = 1:length(ids)
%     filename=files{i};
%     fprintf('File processing: %s\n', filename);
    id = ids(i);
    %%%datalocation = [data_dir '/' num2str(id)];
    datalocation = [data_dir];
    %[posterior, out] = bpd_trust_Qlearning_ushifted(id, counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret);
    [posterior, out] = bpd_trust_Qlearning_ushifted(id, counter, datalocation,...
        'multisession', multisession, 'fixed', fixed_params_across_runs,...
        'sigmakappa', sigma_kappa, 'reputation_sensitive', reputation_sensitive,...
        'humanity', humanity, 'valence_p', valence_p, 'valence_n', valence_n,...
        'assymetry_choices', assymetry_choices, 'regret', regret, ...
        'local_data_dir', local_data_dir);
    L(i) = out.F;
    
    
    b=bpd_trustmakeregressor_group(id);
    
end


L_name = sprintf('L_counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d_regret%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret);
save(char(L_name), 'L'); %Just saveing L's
