function par=trust_extract_PEs_bpd(subdirs)
close all;

%subdirs = {'f_trust_Qlearn_counter_hybrid' 'f_trust_Qlearn_counter_hybrid_regret', 'f_trust_Qlearn_counter_trustee'};
%subdirs = {'f_trust_Qlearn_counter_hybrid'};
%subdirs={'f_trust_Qlearn1', 'f_trust_Qlearn_counter_corrected', 'f_trust_Qlearn_counter_hybrid', 'f_trust_Qlearn_counter_hybrid_regret', 'f_trust_Qlearn_counter_hybrid_regret_purist', 'f_trust_Qlearn_counter_trustee'};

%Really need to fix this cd stuff
current_dir=pwd;


for j = 1:length(subdirs)
    %Quick username check, and path setting, this may have to change depending
    %on the machine you are currently working on!
    os = computer;
    if strcmp(os(1:end-2),'PCWIN')
         datalocation = glob(['E:\data\bpd_trust\vba_out' filesep]);
        %datalocation = glob(['E:\data\hallquist trust' filesep]);
    else
        [~, me] = system('whoami');
        me = strtrim(me);
        if strcmp(me,'polinavanyukov')==1
            datalocation = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/');
        else
            datalocation = glob('?');
        end
    end
    
    %% choose model's parameters 
    par.multisession = 0;
    par.fixed_params_across_runs = 1;
    par.sigma_kappa = 1;
    par.reputation_sensitive = 0;
    par.humanity = 0;
    par.valence_p = 0;
    par.valence_n = 0;
    par.assymetry_choices = 0;
    par.regret = 0;
    
    % get ID list
    cd([datalocation{1} filesep subdirs{j}]);
    
    %get files
    files = dir(strcat('*',sprintf('model*_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d', par.multisession, par.fixed_params_across_runs, par.sigma_kappa, par.reputation_sensitive, par.humanity, par.valence_p, par.valence_n, par.assymetry_choices, par.regret),'.mat'));
    
    %Grab model name -- figure out a away to do this in one line...
    expression  = 'model_(\w+)_m';
    model_name = cellfun(@(x) cell2mat(x),regexp(files(1).name,expression,'tokens'),'UniformOutput',0);
    model_name = model_name{:};
    
    par.([subdirs{j} 'model_name']) = model_name;
    
    num_of_subjects = length(files);
    P=zeros(num_of_subjects,193); %CALCULATED PEs
    M=zeros(num_of_subjects,193); %MODEL PEs
    N=zeros(num_of_subjects,193); %value
    for ct = 1:num_of_subjects
        filename=files(ct).name;
        fprintf('File processing: %s\n', filename);
        %subject_id = filename(isstrprop(filename,'digit'));
        load(filename);
        if ischar(id)
            subject_id = str2double(id);
        else
            subject_id = id;
        end
        P(ct,:) = [subject_id, out.suffStat.PE];
        M(ct,:) = [subject_id, out.suffStat.muX(2,:)];
        N(ct,:) = [subject_id, out.suffStat.muX(1,:)];
        %close all;
    end
    
    if ~exist('PEs', 'dir')
        mkdir('PEs')
    end
    
    
    P_name = ['PEs' filesep sprintf('calcPEs_model_%s_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d',model_name, par.multisession, par.fixed_params_across_runs, par.sigma_kappa, par.reputation_sensitive, par.humanity, par.valence_p, par.valence_n, par.assymetry_choices, par.regret)];
    save(char(P_name), 'P');
    
   % M_name = sprintf('modelPEs_counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d_regret%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret);
    M_name = ['PEs' filesep sprintf('modelPEs_model_%s_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d',model_name, par.multisession, par.fixed_params_across_runs, par.sigma_kappa, par.reputation_sensitive, par.humanity, par.valence_p, par.valence_n, par.assymetry_choices, par.regret)];
    save(char(M_name), 'M');
    
    %N_name = sprintf('values_counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices);
    N_name = ['PEs' filesep sprintf('values_model_%s_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d',model_name, par.multisession, par.fixed_params_across_runs, par.sigma_kappa, par.reputation_sensitive, par.humanity, par.valence_p, par.valence_n, par.assymetry_choices, par.regret)];
    save(char(N_name), 'N');
    close all;
end

cd(current_dir)