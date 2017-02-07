function [posterior,out] = bpd_trust_Qlearning_ushifted(id, datalocation, varargin)
% Model fitting of Qlearning to trust data for BPD study, using VBA-toolbox
% Written by Polina Vanyukov, Alex Dombrovski, & Jonathan Wilson
%
% REQUIRED INPUTS:
% id =           6-digit ID 
% datalocation = Where the subjects .mat file is located
%
% OPTIONAL INPUTS:
% multisession = 1 (default)            treat each trustee as a different run, if 0 all runs are concatenated
% fixed = 1 (default)                   fix learning rate and inv. temperature, BUT NOT X0, across trustees
% reputation_sensitive = 0 (default)    modelling trustees' reputation
% sigmakappa = 0 (default)              kappa (or action bias) parameter
% humanity = 0 (default)                modelling humanity
% valence_p = 0 (default)               modelling valence of trustees
% valence_n = 0 (default)               modelling valence of ...?
% assymetry_choices = 0 (default)       modelling assymetry in choices
% regret = 0 (default)                  modeling regret
% local_data_dir = '' (default)         where to store workspace of script
% model = 'null' (default)              where model to be evaluated (i,e, which reward matrix to use)
%
% OUTPUTS:
% posterior                     posterior distributions
% out                           fit statistics, diagnostics

%Set up input parser 
p = inputParser;
default_multisession = 1;
default_fixed = 1;
default_reputation_sensitive = 0;
default_sigmakappa = 0;
default_humanity = 0;
default_valence_p = 0;
default_valence_n = 0;
default_assymetry_choices = 0;
default_regret = 0;
default_local_data_dir = ''; %NOTE: unless provided workspace will not save
default_graphics = 0;
default_model = 'null'; %Default model will always be NULL

addRequired(p,'id',@isnumeric);
%addRequired(p,'counter',@isnumeric);
addRequired(p,'datalocation',@isstr);
addParameter(p,'multisession',default_multisession,@isnumeric);
addParameter(p,'fixed',default_fixed,@isnumeric);
addParameter(p,'reputation_sensitive',default_reputation_sensitive,@isnumeric);
addParameter(p,'sigmakappa',default_sigmakappa,@isnumeric);
addParameter(p,'humanity',default_humanity,@isnumeric);
addParameter(p,'valence_p',default_valence_p,@isnumeric);
addParameter(p,'valence_n',default_valence_n,@isnumeric);
addParameter(p,'assymetry_choices',default_assymetry_choices,@isnumeric);
addParameter(p,'regret',default_regret,@isnumeric);
addParameter(p,'local_data_dir',default_local_data_dir,@isstr);
addParameter(p,'model',default_model,@isstr);

%Return error when parameters don't match the schema
p.KeepUnmatched = false; 

%Graphics on/off switch
if ~default_graphics
    options.DisplayWin = 0;
    options.GnFigs = 0;
end


%Parse Params
parse(p,id,datalocation,varargin{:})
params = p.Results;

%Clear all open figs
close all
%% Evolution and observation functions

switch params.model
    case 'null'
        f_fname = @f_trust_Qlearn1; % evolution function (Q-learning) with a single hidden state, Q(share)
    case 'subjectCounterfactual'
        f_fname = @f_trust_Qlearn_counter_hybrid;
    case 'subjectCounterfactual_regret'
        f_fname = @f_trust_Qlearn_counter_hybrid_regret;
    case 'trustee'
        f_fname = @f_trust_Qlearn_counter_trustee;
    otherwise
        error('Model not found, see help for more details')
end

if params.assymetry_choices == 1 || params.sigmakappa == 0
    g_fname=@g_trust_softmax1;      % observation function (softmax mapping), evaluates Q(share)
else
    g_fname = @g_trust_softmax_ED;  % observation function (softmax mapping), evaluates Q(share), w/ kappa parameter
end

%2 = track the value of sharing and PE 1 = only track the value of sharing, i.e. V(trustee)
n_hidden_states = 2; 


%ntrials = 192; %Why is number of trials hardcoded?

%% Load subject's data
load([datalocation filesep sprintf('trust%d',id)])

ntrials = length(b.TrialNumber);

%Parse subjects data
if not(exist('decisions'))
    share =~cellfun(@isempty,strfind(b.PartDecides,'share'));
    keep =~cellfun(@isempty,strfind(b.PartDecides,'keep'));
    missed = ~cellfun(@isempty,strfind(b.PartDecides,'noresponse'));
    b.decisions = zeros(ntrials, 1);
    b.decisions(share) = 1;
    b.decisions(keep) = -1;
    b.decisions(missed) = 0;
end
if exist('share')
    actions = share(1:ntrials)'; %subject's actions
else
    share = zeros(length(b.decisions),1);
    share(b.decisions==1) = 1;
    actions = share(1:ntrials);
end
if not(exist('noresponse'))
    noresponse = zeros(length(b.decisions),1);
    noresponse(b.decisions==0|b.decisions==-999)=1;
end
actions = double(actions);
u(1,:) = actions;

rewards = double(strcmp(b.TrusteeDecides,'share')); %rewards including counterfactual ones (trustee's actions)
rewards(rewards==0) = -1;
u(2,:) = rewards(1:ntrials)';

%% Experimental design
%reputation vector
trustee_ID = zeros(length(b.identity),1);
trustee_ID(strcmp(b.identity,'good')) = 1;  % Lina's initial coding = 2
trustee_ID(strcmp(b.identity,'bad')) = -1; 
trustee_ID(strcmp(b.identity,'neutral')) = 0; % Lina: -1
u(3,:) =  trustee_ID(1:ntrials);

% humanity vector
human = ones(length(b.identity),1);
human(strcmp(b.identity,'computer')) = 0;
u(5,:) = human(1:ntrials);

% positive/negative valence vector
valence_pos = zeros(length(b.identity),1);
valence_neg = zeros(length(b.identity),1);
valence_pos(strcmp(b.identity,'good')) = 1;
valence_neg(strcmp(b.identity,'bad')) = 1;
u(6,:) = valence_pos(1:ntrials);
u(7,:) = valence_neg(1:ntrials);

%initial state value only to be sensitive
index_vector = zeros(length(b.identity),1);
blocklength = 48;
index_vector([1,1+blocklength,1+2*blocklength, 1+3*blocklength],1) = 1;
u(4,:) = index_vector(1:ntrials);

y = u(1,:); %the subject's actions
%shifting u
u = [zeros(7,1) u(:,1:end-1)];


%% Sensitivities or bias parameters
options.inF.reputation_sensitive = params.reputation_sensitive;
options.inF.humanity = params.humanity;
options.inF.valence_p = params.valence_p;
options.inF.valence_n = params.valence_n;
options.inF.assymetry_choices = params.assymetry_choices;
options.inF.regret = params.regret;

%% dimensions and options
n_t = size(u(2,:),2); % number of trials
n_trustees = 4;
if params.multisession
    options.multisession.split = repmat(n_t/n_trustees,1,n_trustees); % two sessions of 120 datapoints each
    %% fix parameters
    if params.fixed
        options.multisession.fixed.theta = 2; %fixing sensitivity parameter to be the same across sessions
        options.multisession.fixed.phi = 'all';
        options.multisession.fixed.X0 = 'all';
    end
        
end
options.isYout=zeros(size(u(1,:)));
options.isYout(:,1)=1;
options.isYout(noresponse'==1)=1;
options.binomial = 1;
options.skipf = zeros(1,n_t);
options.skipf(1) = 1; % apply identity mapping from x0 to x1.

%% defined number of hidden states and parameters
n_theta = 1+params.reputation_sensitive+params.humanity+params.valence_p+params.valence_n+params.assymetry_choices+(params.regret); %evolution function paramaters: by default = 1 (learning rate), adds 1 for each additional experimental design element
n_phi = 1+params.sigmakappa; %observation function parameters: by default = 1 (beta), adds 1 for each additional observation parameter
dim = struct('n',n_hidden_states,'n_theta',n_theta,'n_phi',n_phi, 'n_t', n_t);


%% priors
if n_phi ==2
    priors.muPhi = [1;0];
else 
    priors.muPhi = 1;
end

priors.muTheta = zeros(dim.n_theta,1);
priors.muX0 = zeros(n_hidden_states,1);
if params.reputation_sensitive||params.humanity||params.valence_p||params.valence_n
%    priors.SigmaTheta = 1e1*eye(dim.n_theta);     
    priors.SigmaTheta = diag([10 1/3*ones(1,dim.n_theta-1)]);
elseif params.assymetry_choices 
    priors.SigmaTheta = diag([10 10]);
elseif params.regret
    priors.SigmaTheta = diag([10 1/3*ones(1,dim.n_theta-1)]);
else
    priors.SigmaTheta = 10;
end

if params.sigmakappa
    priors.SigmaPhi = diag([10 params.sigmakappa/3]);
else
    priors.SigmaPhi = 10;
end

%% set priors on initial states
%priors.SigmaX0 = .3*eye(dim.n);% tracking a single hidden state (value)
priors.SigmaX0 = diag([.3 0]);  % tracking value and prediction error
%priors.SigmaX0 = 0*eye(dim.n); %used this before for tacking value and
%prediction error
% priors.SigmaX0 = zeros(dim.n); % because of correlation with learning
% % rates, one may want to fix initial states to analyze trustee-wise
% % learning rates

%% set hyper-priors
% priors.a_sigma = 1;       % Jeffrey's prior
% priors.b_sigma = 1;       % Jeffrey's prior
% priors.a_alpha = Inf;
% priors.b_alpha = 0;

options.priors = priors;

options.verbose=1;
%options.DisplayWin=0;
options.GnFigs=0;

%% model inversion
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
%displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf);

%% print condition order for interpreting results
ConditionOrder  = unique(b.identity,'stable');
out.design = ConditionOrder;

%h = figure(1);
%savefig(h,sprintf('%d_counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choices%d_regret%d', id, counter, multisession, fixed, sigmakappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret));

%% get prediction errors
alpha = 1./(1+exp(-posterior.muTheta(1)));
if ~params.multisession
    out.suffStat.PE = diff(out.suffStat.muX)./alpha;
else
    out.suffStat.PE = diff(sum(out.suffStat.muX))./alpha;
end

if not(isnumeric(id))
    id = str2double(id);
end

%Only save is local_data_dir param is provided
if ~isempty(params.local_data_dir)
    %Ensure proper syntax
    if ~(strcmp(params.local_data_dir(end),'/') || strcmp(params.local_data_dir(end),'\'))
        params.local_data_dir = [params.local_data_dir filesep];
    end
    
    subdir = func2str(f_fname);
    
    if ~exist([params.local_data_dir subdir],'dir')
        mkdir([params.local_data_dir subdir])
    end
    %Save all params in file name
    filename = sprintf('bpd_trust_%s_cntr%d_mltrun%d_fixed%d_kappa%d_rep%d_hum%d_val_p%d_val_n%d_as_choices%d_reg%d', id,...
        params.model,params.multisession, params.fixed, params.sigmakappa, params.reputation_sensitive, params.humanity,...
        params.valence_p, params.valence_n, params.assymetry_choices, params.regret);
    save([params.local_data_dir filename]);
end






 