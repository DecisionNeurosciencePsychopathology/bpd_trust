function b = bpd_trustmakeregressor_group(id)
%Grabbing all subjects that have already been processed and turned 
%into .mat file. Then making regressors for each one.

%grabbing files
%data_dir_str= '/Volumes/bek/trust_analyses/regs/scan_behavior/';
%data_dump_str = '/Volumes/bek/trust_analyses/regs/';
%grabbing files
data_dir_str= sprintf('subjects/%d',id);
data_dump_str = sprintf('regs/%d/',id);

if ~exist(data_dump_str,'dir')
    mkdir(data_dump_str)
    fprintf('Creating id specific reg folder in: %s\n\n',data_dump_str);
end

%cd(data_dir_str)
files = glob([data_dir_str sprintf('/trust*.mat',id)]);


filename = char(files);
num_of_subjects = 1;

%hard coded times
dur_choice_display = 300;
dur_feedback = 1200; %1100msec in the file?

%scanning parameters
scan_tr = .6;

%deteremine number of blocks/trials per block
rsf=dir([data_dir_str sprintf('/*_scan_*.edat2',id)]);
rsfname=rsf.name;
run_asterisk =0; %Default moving forward


if strfind(rsfname,'072916')>=1
    trialsperblock= 192;
    blocknum= 1;
    block_length= 1685;
    
else
    run_asterisk =1;
    trialsperblock= 48;
    blocknum= 4;
    block_length= 440;
    
end

%Might need if statement here...
%block_length = 155; %This will depend on subject if before 7/26

block_end = block_length*scan_tr*1000; %in msec
hemoir = spm_hrf(scan_tr, [6,16,1,1,6,0,32]); % better than resampling and smoothing 
frequency_scale_hz = 10;
% this scale is in msec, but it is separated into bins of X
% Hz (defined by 'frequency_scale' above. the resulting
% output will be in the scale of X Hz.
bin_size = 1/frequency_scale_hz*1000; % convert Hz to ms

%for index = 44
%for index=10
for index=1:num_of_subjects
    %filename = files(index).name;
    fprintf('File processing: %s\n', filename);
    %id = filename(isstrprop(filename,'digit'));
    if not(str2double(id)==219956||str2double(id)==220017) %Remove later!
        load(filename);
        data_dump_str = strcat(sprintf('regs/%d/',id), num2str(id));
        
        if run_asterisk==1
           data_dump_str= strcat(sprintf('regs/%d/',id), num2str(id));
        end
        
        b.regs = [];    
        
    %Subjects may be using different trial versions of the task, 
    %which only affects how the timings for the scanner should be calculated. 
    %For some participants first fixation onset is the ITI fixation(1). 
    %For others, first fixation onset is labeled as such, 2nd fixation onset 
    %= ITI fixation (1), etc.
        if not(exist('decisions'))
            share =~cellfun(@isempty,strfind(b.PartDecides,'share'));
            keep =~cellfun(@isempty,strfind(b.PartDecides,'keep'));
            noresponse = ~cellfun(@isempty,strfind(b.PartDecides,'noresponse'));
            b.decisions = zeros(192, 1);
            b.decisions(share) = 1;
            b.decisions(keep) = -1;
            b.decisions(noresponse) = 0;
        end
        
        if length(b.ITIfixation_OnsetTime) < 192
            b.ITIfixation_OnsetTime(length(b.ITIfixation_OnsetTime)+1:192)=-999;
            b.ITIfixation_OffsetTime(145:192)=num2cell(-999);
        end
        
        if b.ITIfixation_OnsetTime(192) == -999 || not(iscell(b.firstFixation_OnsetTime))
            firstfix_Onset = b.firstFixation_OnsetTime(1);
            trial2_ITI = 1;
            trial48_ITI = trialsperblock-1;
        else firstfix_Onset = b.ITIfixation_OnsetTime(1);
            trial2_ITI = 2;
            trial48_ITI = trialsperblock;
        end

        %for those participants for whom the partnerchoice offset time was not
        %recorded by e-prime
        if iscell(b.partnerchoice_OffsetTime)
            partnerchoice_OffsetTime=b.displaychoice_OnsetTime-1;
        else
            partnerchoice_OffsetTime=b.partnerchoice_OffsetTime;
        end

        fixations = [];
        fixations.event_beg=zeros(trialsperblock,blocknum);
        fixations.event_end=zeros(trialsperblock,blocknum); 
        taskness.event_beg =zeros(trialsperblock,blocknum);
        taskness.event_end=zeros(trialsperblock,blocknum);
        decision.event_beg=zeros(trialsperblock,blocknum);
        decision.event_end=zeros(trialsperblock,blocknum);
        feedback.event_beg=zeros(trialsperblock,blocknum);
        feedback.event_end=zeros(trialsperblock,blocknum);
        b.missed_trials = (b.partnerchoice_RESP==-999);
        b.notmissed_trials = (b.missed_trials-1)*-1; 
        %plot(b.notmissed_trials);
        trial1_index = 1;
        trial_index = trialsperblock;
        
        for block= 1:blocknum
            %for fixation screens
            if block == 1 && (b.ITIfixation_OnsetTime(192) == -999 || not(isempty([b.firstFixation_OnsetTime])))
                fixations.event_beg(:,block) = [0; b.ITIfixation_OnsetTime(trial2_ITI:trial48_ITI)-firstfix_Onset];
            else
                %fixations.event_beg(:,block) = b.ITIfixation_OnsetTime(trial2_ITI-1:trial48_ITI)-firstfix_Onset+block_end*(block-1);
                fixations.event_beg(:,block) = b.ITIfixation_OnsetTime(trial2_ITI-1:trial48_ITI)-firstfix_Onset;
            end        
            %fixations.event_end(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial48_index) - firstfix_Onset+block_end*(block-1);
            fixations.event_end(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial_index) - firstfix_Onset;

            %for trial onset to offset; Taskness
            taskness.event_beg(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial_index)-firstfix_Onset;
            taskness.event_end(:,block) = b.outcome_OffsetTime(trial1_index:trial_index)-firstfix_Onset;

            %for decision onset to response (motor response)
            decision.event_beg(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial_index)-firstfix_Onset;
            decision.event_end(:,block) = partnerchoice_OffsetTime(trial1_index:trial_index)-firstfix_Onset; 

            %for feedback onset to offset
            feedback.event_beg(:,block) = b.outcome_OnsetTime(trial1_index:trial_index)-firstfix_Onset;
            feedback.event_end(:,block) = b.outcome_OffsetTime(trial1_index:trial_index)-firstfix_Onset; 
            
            %% New for aligning value: RT aligned to feedback
%             decision.event_beg(:,block) = b.partnerchoice_RTTime(trial1_index:trial48_index)-firstfix_Onset;
%             decision.event_end(:,block) = partnerchoice_OffsetTime(trial1_index:trial48_index)-firstfix_Onset; 

            %epoch window + missed trials + to censor regressor 
            epoch_window = 0:bin_size:taskness.event_end(48, block);
            %% There is a strange glitch, likely owing to the createSimpleRegressor code;
            %% modified by AD to include ITI
            if(sum(b.missed_trials(trial1_index:trial_index))) > 0
                % write 1s for missed trials to censor
                tmp_reg.(['regressors' num2str(block)]).to_censor = createSimpleRegressor(taskness.event_beg,taskness.event_end, epoch_window, b.missed_trials(trial1_index:trial_index));
            else
                % write a vector of 0s the size of regressors
                tmp_reg.(['regressors' num2str(block)]).to_censor = zeros(size(createSimpleRegressor(taskness.event_beg,taskness.event_end, epoch_window, b.notmissed_trials(trial1_index:trial_index))));
            end
            % flip to_censor
            tmp_reg.(['regressors' num2str(block)]).to_censor = 1-tmp_reg.(['regressors' num2str(block)]).to_censor;

            %% GSR resample
            tmp = gsresample( ...
             [zeros(50,1)' tmp_reg.(['regressors' num2str(block)]).to_censor(1:end-51)], ...
             10,1./scan_tr);
            tmp = ceil(tmp);
            tmp = [tmp ones(1, (block_length-1)-length(tmp))];
            tmp = [tmp zeros(1,155-length(tmp))];
            tmp_reg.(['regressors' num2str(block)]).to_censor = tmp;
            %plot(tmp_reg.(['regressors' num2str(block)]).to_censor);
            
            %% Censoring blocks w/ movement
            tmp_reg=censorMovement(id, tmp_reg, block);            
            %plot(tmp_reg.(['regressors' num2str(block)]).to_censor);
            
            %for next loop iteration, reinitilize variables
            if block < 4
                firstfix_Onset = b.ITIfixation_OnsetTime(trial2_ITI-1+48);
            end
            trial2_ITI=trial2_ITI+trialsperblock;
            trial48_ITI=trial48_ITI+trialsperblock;
            trial1_index = trial1_index+trialsperblock;
            trial_index = trial_index+trialsperblock;       
        end
        fixations.event_beg=reshape(fixations.event_beg,[192,1]);
        fixations.event_end=reshape(fixations.event_end,[192,1]);
        taskness.event_beg =reshape(taskness.event_beg,[192,1]);
        taskness.event_end=reshape(taskness.event_end,[192,1]);
        decision.event_beg=reshape(decision.event_beg,[192,1]);
        decision.event_end=reshape(decision.event_end,[192,1]);
        feedback.event_beg=reshape(feedback.event_beg,[192,1]);
        feedback.event_end=reshape(feedback.event_end,[192,1]);
        %plotting to censor regressors
%         plot(tmp_reg.regressors3.to_censor);
%         figure(2);
%         plot(b.notmissed_trials(192-48:192));
        %concatenating
        if blocknum>1
            b.to_censor = [tmp_reg.regressors1.to_censor tmp_reg.regressors2.to_censor tmp_reg.regressors3.to_censor tmp_reg.regressors4.to_censor]; 
        else
                b.to_censor=tmp_reg.regressors1.to_censor;
        end
 %       plot(b.to_censor);
        b.to_censor = transpose(b.to_censor);
        
        %exclude missed trials from other regressors
        b.notmissed_trials = (b.partnerchoice_RESP~=-999);
%           [b.stim_times.fix_fsl,b.stim_times.fix_spmg]=write3Ddeconv_startTimes(data_dump_str,fixations.event_beg,fixations.event_end,'fixation_Times',b.notmissed_trials,0);
%           [b.stim_times.task_fsl,b.stim_times.task_spmg]=write3Ddeconv_startTimes(data_dump_str,taskness.event_beg,taskness.event_end,'task_Times',b.notmissed_trials,0);
%           [b.stim_times.resp_fsl,b.stim_times.resp_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,decision.event_end,'decision_Times',b.notmissed_trials,0);
  
        % right = 2; left = 7;
        b.leftVSright = zeros(size(b.partnerchoice_RESP));
        b.leftVSright(b.partnerchoice_RESP==7)=1;
        b.leftVSright(b.partnerchoice_RESP==2)=-1; 

        % subject's decisions, share/keep;
        b.shareVSkeep = zeros(size(b.partnerchoice_RESP));
        b.shareVSkeep(b.decisions==1 & b.partnerchoice_RESP ~= -999) = 1;
%        b.shareVSkeep(b.decisions==-1 & b.partnerchoice_RESP ~= -999) = -1;

%           [b.stim_times.left_fsl,b.stim_times.left_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,decision.event_end,'leftVSright',b.leftVSright,0);
%           [b.stim_times.share_fsl,b.stim_times.share_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,decision.event_end,'p_shareVSkeep',b.shareVSkeep,0);
           [b.stim_times.share_fsl,b.stim_times.share_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'p_shareVSkeep',b.shareVSkeep,0);
      
        share =~cellfun(@isempty,strfind(b.TrusteeDecides(1:192),'share'));
        keep =~cellfun(@isempty,strfind(b.TrusteeDecides(1:192),'keep'));
        b.t_share = zeros(192,1);
        b.t_keep = zeros(192,1);
        b.t_share(share) = 1;
        b.t_share(keep) = -1;

%           [b.stim_times.feedback_fsl,b.stim_times.feedback_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'feedback_Times',b.notmissed_trials,0);
%           [b.stim_times.t_shared_fsl,b.stim_times.t_shared_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'tshared_Times',b.t_share,0);

       
        %% trustee contrasts 
        b.trustee_GB = zeros(length(b.identity),1);
        b.trustee_G0 = zeros(length(b.identity),1);
        b.trustee_HC = ones(length(b.identity),1);
%        b.trustee_BA = ones(length(b.identity),1)*-1;
        b.trustee_BA = zeros(length(b.identity),1);
        b.trustee_survey_GB = zeros(length(b.identity),1);
        b.trustee_action = ones(length(b.identity),1);
                
%         %survey elicited scaling (pre-rating only)
%         load('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/surveys/trust_group_survey.mat');
%         if ischar(id)
%             rows = group_survey_table.ID==str2double(id);
%         else
%             rows = group_survey_table.ID==id;
%         end
%         survey = group_survey_table(rows,:);
%         good_rating = survey.Rating_Pre(strcmp(survey.Trustee,'good') & strcmp(survey.Question,'trustworthy'));
%         bad_rating = survey.Rating_Pre(strcmp(survey.Trustee,'bad') & strcmp(survey.Question,'trustworthy'));
%         %neutral_rating = group_survey_table.Rating_Pre(group_survey_table.ID == id && strcmp(group_survey_table.Trustee,'neutral') && strcmp(group_survey_table.Question,'trustworthy'));
%         %computer_rating = group_survey_table.Rating_Pre(group_survey_table.ID == id && strcmp(group_survey_table.Trustee,'computer') && strcmp(group_survey_table.Question,'trustworthy'));
%         b.trustee_survey_GB(strcmp(b.identity,'good')) = 1*good_rating/7; 
%         b.trustee_survey_GB(strcmp(b.identity,'bad')) = -1*bad_rating/7;
%        [b.stim_times.surveyGB_fsl,b.stim_times.surveyGB_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'surveyGB',b.trustee_survey_GB,0);
        
        % good VS bad
        b.trustee_GB(strcmp(b.identity,'good')) = 1;         
        b.trustee_G0(strcmp(b.identity,'good')) = 1;
        b.trustee_GB(strcmp(b.identity,'bad')) = -1;
        
        [b.stim_times.trusteeGB_fsl,b.stim_times.trusteeGB_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'trusteeGB',b.trustee_GB,0);        
        [b.stim_times.trusteeG0_fsl,b.stim_times.trusteeG0_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'trusteeG0',b.trustee_G0,0);        

        % human vs non-human
        b.trustee_HC(strcmp(b.identity,'computer')) = 0;
       [b.stim_times.trusteeHC_fsl,b.stim_times.trusteeHC_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'trusteeHC',b.trustee_HC,0);
%         
       % bad vs rest
       b.trustee_BA(strcmp(b.identity,'bad')) = 1;
%       [b.stim_times.trusteeBA_fsl,b.stim_times.trusteeBA_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'trusteeBA',b.trustee_BA,0);        
            
       % trustee action regressor
       b.trustee_action(strcmp(b.TrusteeDecides,'keep'))=-1;
%       [b.stim_times.trustee_action_fsl,b.stim_times.trustee_action_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'trusteeAction',b.trustee_action,0);
       
       
       % GB contrast needs a censor regressor that censors other trustee
       % blocks
       b.trusteeBYtactionGB = b.trustee_GB.*b.trustee_action;
%       [b.stim_times.trusteeBYactionGB_fsl,b.stim_times.trusteeBYactionGB_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'GBxtaction',b.trusteeBYtactionGB,0);
       b.trusteeBAxp_shareVSkeep = b.trustee_BA.*b.shareVSkeep;
%       [b.stim_times.trusteeBYactionBA_fsl,b.stim_times.trusteeBYactionBA_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'BAxp_decision',b.trusteeBAxp_shareVSkeep,0);
       
       b.trusteeGBxp_shareVSkeep = b.trustee_GB.*b.shareVSkeep;
       [b.stim_times.trusteeBYactionGB_fsl,b.stim_times.trusteeBYactionGB_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'GBxp_decision',b.trusteeGBxp_shareVSkeep,0);
       
       b.trusteeG0xp_shareVSkeep = b.trustee_G0.*b.shareVSkeep;
       [b.stim_times.trusteeBYactionG0_fsl,b.stim_times.trusteeBYactionG0_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'G0xp_decision',b.trusteeG0xp_shareVSkeep,0);
  
       b.trusteeHCxp_shareVSkeep = b.trustee_HC.*b.shareVSkeep;
       [b.stim_times.trusteeBYactionHC_fsl,b.stim_times.trusteeBYactionHC_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'HCxp_decision',b.trusteeHCxp_shareVSkeep,0);
 
%       b.surveyGBbytactionGB = b.trustee_survey_GB * b.trustee_action;
%         gdlmwrite(strcat(data_dump_str, 'to_censor'),[b.to_censor],'\t');
         gdlmwrite(strcat(data_dump_str, 'to_censor_motion_corr'),[b.to_censor],'\t');          
%         gdlmwrite(strcat(data_dump_str, 'to_censor_PEs'),[b.to_censor_PEmaps],'\t');
    end
end

%call asterisk if there were 4 blocks
if run_asterisk== 1
    
    asterisk(id);
    
    rmdir(sprintf('%d',id),'s')
    
    cd ..
    
end



return


