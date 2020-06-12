%% OVERVIEW
% ===============================

% The cardiac model was a simple attempt at just using ellipses
% Motion is modelled by sinusoid behaviour for both cardiac and respiratory
% Cardiac motion is a "one-directional" sine (just -ve, 0 for +ve)
% The tissue parameters are coded in to the simulator

% ===============================
%
% Requirements and defaults can be found by
help sim_SL_cardiac;

% I dont recommend running this script! It includes examples to show some configurations..

%% REQUIREMENTS
% ===============================

% Note that this simulate requires

% ------------------------------
% MRIPhantomv0-8
% ------------------------------
    % from Matthieu Guerquin-Kern @ Lausanne, France
    % https://miplab.epfl.ch/pub/guerquinkern1001.pdf
%
% available at
% http://bigwww.epfl.ch/algorithms/mriphantom/
%
% and to make the mex files (for coil estimation)
% add all paths:
    % mriphantom_path = uigetdir(matlabroot, 'Select sim dir');
    % addpath(genpath([mriphantom_path filesep 'MRIPhantomv0-8']))

% ------------------------------
% The ISMRM sunrise toolbox
% ------------------------------
    % from Michael Hansen @ NIH, USA
    % http://hansenms.github.io/sunrise/
%
% available at
% https://github.com/hansenms/ismrm_sunrise_matlab
%
% with the following paths added:
    % ismrm_sunrise_path = uigetdir(matlabroot, 'Select ISMRM sunrise dir');
    %
    %     addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master'])
    %     addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'vdspiral'])
    %     addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'utilities'])
    %     addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'graph'])
    %     addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'systems'])
    %     addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'wls'])
    %     addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'penalty'])
    %     addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'general'])
    %     addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'fbp'])
    %     addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'nufft'])
    %     addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'mri'])
    %     addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'mex' filesep 'v7']) 
    %     addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'mri'])

% ===============================

%% DEFAULTS
% ===============================

% Choose from:  
%   
% sampling design 
% enc_info.res            = [256 256]; 
% enc_info.FOV            = [32 32]; 
% enc_info.k              = GenerateFullCart2DKspace(enc_info.res,enc_info.FOV); % Arbitrary trajectory (everything nuffted) 
% enc_info.k              = reshape(permute(enc_info.k,[2 1]),[enc_info.res 2]); % nufft format : [samples shots 2] 
% enc_info.coils          = 6; 
% enc_info.rel_coil_snr   = [];                          % defined relative coil SNR 
% enc_info.sens_maps      = []; 
%   
% contrast 
% signal_info.TR              = 10; 
% signal_info.numAcqs         = 24; 
% signal_info.contrast        = 1; 
% signal_info.time_var_sig    = 0; 
% signal_info.inversion_alpha = 1.8; 
% signal_info.image_snr       = 50; 
% signal_info.snr_scale       = [];                      % custom/consistent scale factor to transform 
% signal_info.ext_sig_flag    = 0;                       % [0]: simple model, [1]: arbitrary signal model 
% signal_info.ext_sig_profs   = [];                      % [NumAcqs x 5] : [fat muscle blood liver scar] 
%   
% motion 
% model_motion.cardiac_interval   = 777;                 % ms 
% model_motion.cardiac_variation  = 6;                   % percentage 
% model_motion.cardiac_phase      = randn(1);            % offset to starting phase 
% model_motion.resp_interval      = 3318;                % ms 
% model_motion.resp_variation     = 13;                  % percentage 
% model_motion.resp_phase         = randn(1);            % offset to starting phase 
% model_motion.motion_model       = @(t_ms, int_ms, phase) sin(2*pi*((t_ms/int_ms) + phase)); 
  

% ===============================

% ===============================
%% EXAMPLES
% ===============================

% =============
% quick test (Cartesian)
% =============
enc_info.coils           = 0; % homogeneous coil
signal_info.numAcqs      = 10;
signal_info.image_snr    = 100; % first run came out really low!
signal_info.time_var_sig = 1; % run inversion prep (just one for now)
model_motion.dummy       = 1; % just dummy interface parameter for the moment

[kdata, sim_info] = sim_SL_cardiac(enc_info, signal_info, model_motion);

% =============
%% quick test (spiral)
% =============
    interleaves = 24; matrix = 256; FOV = 32; VDSf = 100;
    gmax = 2.4; smax = 14414.4; dt = 2e-6; krmax = 1/(2*(FOV(1)/matrix ));
    [k,g] = vds(smax,gmax,dt,interleaves,FOV,krmax); close;
    samples = length(k); trajectory_nominal = zeros(samples,interleaves,2);

    for i = 1:interleaves
        rot = 2*pi*((double(i)-1)/double(interleaves));
        trajectory_nominal(:,i,1) = -( real(k(1:samples)) *cos(rot) + imag(k(1:samples)) *sin(rot));
        trajectory_nominal(:,i,2) = -(-real(k(1:samples)) *sin(rot) + imag(k(1:samples)) *cos(rot));
    end

enc_info.k               = trajectory_nominal; % homogeneous coil

enc_info.coils           = 0; % homogeneous coil
signal_info.numAcqs      = 30;
signal_info.image_snr    = 100; % first run came out really low!
signal_info.time_var_sig = 1; % run inversion prep (just one for now)
model_motion.dummy       = 1; % just dummy interface parameter for the moment

[kdata, sim_info] = sim_SL_cardiac(enc_info, signal_info, model_motion);

% =============
%% quick test (spiral) with external signal
% =============

signal_info.numAcqs      = 300;
signal_info.image_snr    = 100; % first run came out really low!
signal_info.time_var_sig = 1; % run inversion prep (just one for now)
signal_info.ext_sig_flag = 1;
signal_info.ext_sig_profs = [linspace(0.4,0.8,300)' .25+0.025*abs(randn(1,300)') linspace(0.9,0.5,300)' 0.5*ones(1,300)' .5+0.15*abs(randn(1,300)')];
model_motion.dummy       = 1; % just dummy interface parameter for the moment

[kdata, sim_info] = sim_SL_cardiac(enc_info, signal_info, model_motion);


% =============
%% multi-coil test (sloooow)
% =============
clear enc_info signal_info
enc_info.dummy           = 1; % just dummy parameter for the moment
signal_info.numAcqs      = 10;
signal_info.image_snr    = 100; % first run came out really low!
signal_info.time_var_sig = 0; % run inversion prep (just one for now)
model_motion.dummy       = 1; % just dummy interface parameter for the moment

[kdata, sim_info] = sim_SL_cardiac(enc_info, signal_info, model_motion);

% =============
%% run without motion:
% =============
clear model_motion
model_motion.cardiac_interval = 0;
model_motion.resp_interval = 0;

[kdata, sim_info] = sim_SL_cardiac(enc_info, signal_info, model_motion);

