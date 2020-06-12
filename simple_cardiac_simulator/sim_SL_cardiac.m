function [kdata, sim_info] = sim_SL_cardiac(enc_info, signal_info, model_motion)
% function [kdata, sim_info] = sim_SL_cardiac(enc_info, signal_info, model_motion)
% <a href="matlab:sim_SL_cardiac">run</a> for help description

if nargin < 1
    input_help; % print out input options
%     if nargout > 0
%     end
    return;
end

%% User options with defaults
if exist('enc_info','var')
    enc_info_in = enc_info;
else
    enc_info_in = struct();
end
enc_info_in_list = fieldnames(enc_info_in);

if exist('signal_info','var')
    signal_info_in = signal_info;
else
    signal_info_in = struct();
end
signal_info_in_list = fieldnames(signal_info_in);

if exist('model_motion','var')
    model_motion_in = model_motion;
else
    model_motion_in = struct();
end
model_motion_in_list = fieldnames(model_motion_in);

%% defaults
% sampling design
enc_info.res            = [256 256];
enc_info.FOV            = [32 32];
enc_info.k              = GenerateFullCart2DKspace(enc_info.res,enc_info.FOV);
enc_info.k              = reshape(permute(enc_info.k,[2 1]),[enc_info.res 2]); % nufft format : [samples shots 2]
enc_info.coils          = 6;
enc_info.rel_coil_snr   = [];
enc_info.sens_maps      = [];

% contrast 
signal_info.TR              = 10;
signal_info.numAcqs         = 24;
signal_info.contrast        = 1;
signal_info.time_var_sig    = 0;
signal_info.inversion_alpha = 1.8;
signal_info.image_snr       = 50;
signal_info.snr_scale       = [];
signal_info.ext_sig_flag    = 0;
signal_info.ext_sig_profs   = [];

% motion
model_motion.cardiac_interval   = 777;      % ms
model_motion.cardiac_variation  = 6;       % percentage
model_motion.cardiac_phase      = randn(1); % offset to starting phase
model_motion.resp_interval      = 3318;     % ms
model_motion.resp_variation     = 13;       % percentage
model_motion.resp_phase         = randn(1); % offset to starting phase
model_motion.motion_model       = @(t_ms, int_ms, phase) sin(2*pi*((t_ms/int_ms) + phase));

%% overwrite with user parameters

for i = 1:length(enc_info_in_list)
    if isfield(enc_info, enc_info_in_list{i})
        enc_info.(matlab.lang.makeValidName(enc_info_in_list{i})) = enc_info_in.(matlab.lang.makeValidName(enc_info_in_list{i}));
    end
end


for i = 1:length(signal_info_in_list)
    if isfield(signal_info, signal_info_in_list{i})
        signal_info.(matlab.lang.makeValidName(signal_info_in_list{i})) = signal_info_in.(matlab.lang.makeValidName(signal_info_in_list{i}));
    end
end


for i = 1:length(model_motion_in_list)
    if isfield(model_motion, model_motion_in_list{i})
        model_motion.(matlab.lang.makeValidName(model_motion_in_list{i})) = model_motion_in.(matlab.lang.makeValidName(model_motion_in_list{i}));
    end
end

disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');
disp('<strong>Simulation options:</strong>');                 disp(' ');
disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');
disp('<strong>Encoding:</strong>');                           disp(' ');
disp(enc_info);
disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');
disp('<strong>Signal:</strong>');                           disp(' ');
disp(signal_info);
disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');
disp('<strong>Motion:</strong>');                           disp(' ');
disp(model_motion);
disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');

clear enc_info_in       enc_info_in_list ...
      signal_info_in    signal_info_in_list ...
      model_motion_in   model_motion_in_list;

%% important stuff
k = enc_info.k;
[samples, shots, ~] = size(k);

% v1
t_signal = [0:(signal_info.numAcqs-1)]*signal_info.TR;

% v2
% t_signal = [0:(signal_info.numAcqs-1)]*shots*signal_info.TR;

%% tissue parameters (internal signal profiles)
% overhead for external signal call
% These values are made up for demo purposes. 
tissue_params.fat.M0  = .9; 
tissue_params.fat.T1  = 150;
tissue_params.fat.T2  = 40;
tissue_params.fat.T2s = 20;

tissue_params.muscle.M0  = .85;%0.4; 
tissue_params.muscle.T1  = 1100;
tissue_params.muscle.T2  = 40;
tissue_params.muscle.T2s = 20;

tissue_params.blood.M0  = 1;%0.75; 
tissue_params.blood.T1  = 500;
tissue_params.blood.T2  = 70;
tissue_params.blood.T2s = 40;

tissue_params.liver.M0  = 0.9;%0.3; 
tissue_params.liver.T1  = 800;
tissue_params.liver.T2  = 30;
tissue_params.liver.T2s  = 15;

tissue_params.scar.M0  = 0.9; 
tissue_params.scar.T1  = 400;
tissue_params.scar.T2  = 30;
tissue_params.scar.T2s = 20;

if signal_info.ext_sig_flag == 0
% constant signals:
if signal_info.contrast        == 1
%           % TE | TR | FA | scan
% bssfp
scan_params = [1 signal_info.TR 70 1]; % this should be fed from user
% steady state signal (Bernstein Handbook)
ss_signal = @(tissue_s, scan_params) tissue_s.M0.*( sind(scan_params(3))*(1 - exp(-scan_params(2)./tissue_s.T1)) ).*exp(-scan_params(1)./tissue_s.T2) ...
                                   ./(1 - cosd(scan_params(3))*(exp(-scan_params(2)./tissue_s.T1) - exp(-scan_params(2)./tissue_s.T2)) ... 
                                                            - (exp(-scan_params(2)./tissue_s.T1).*exp(-scan_params(2)./tissue_s.T2)) );
elseif signal_info.contrast      == 2
%           % TE | TR | FA | scan
% GRE
scan_params = [1 signal_info.TR 12 2]; % this should be fed from user
% steady state signal (Bernstein Handbook)
ss_signal = @(tissue_s, scan_params) tissue_s.M0.*( sind(scan_params(3))*(1 - exp(-scan_params(2)./tissue_s.T1)) ).*exp(-scan_params(1)./tissue_s.T2s) ...
                                        ./(1 - cosd(scan_params(3))*exp(-scan_params(2)./tissue_s.T1));
end
                                                        
% time varying signals: (PROBABLY NEED TO EXTERNALISE THIS)
if signal_info.time_var_sig    == 0
    tv_signal = @(tissue_s, scan_params, t_signal, inversion_alpha) abs(ss_signal(tissue_s,scan_params))*ones(size(t_signal));
elseif signal_info.time_var_sig    == 1  
    tv_signal = @(tissue_s, scan_params, t_signal, inversion_alpha) abs(ss_signal(tissue_s,scan_params).*(1 - inversion_alpha*exp(-t_signal/tissue_s.T1)));

% temp = [tv_signal(tissue_params.muscle, scan_params, t_signal, inversion_alpha)' tv_signal(tissue_params.scar, scan_params, t_signal, inversion_alpha)'];
% figure, plot(temp)
end

signal_info.ss_signal   = ss_signal;
signal_info.tv_signal   = tv_signal;
signal_info.scan_params = scan_params;
end

%% Motion vec

% cardiac
current_time = 0;
cMo = zeros(size(t_signal));
if model_motion.cardiac_interval == 0
    % Leave cMo as zeros, no cardiac motion
else
    % Simple half-sinusoid cardiac motion, with cardiac variation
    while (current_time < t_signal(end))
        
        time_index_1 = find(t_signal <= current_time, 1, 'last');
        
        temp_interval = round(model_motion.cardiac_interval + randn(1)*(model_motion.cardiac_variation*model_motion.cardiac_interval/100));
        current_time = current_time + temp_interval;
        
        time_index_2 = find(t_signal <= current_time, 1, 'last');
        t_signal_temp = t_signal(time_index_1:time_index_2) - t_signal(time_index_1);
        temp_interval_motion = model_motion.motion_model(t_signal_temp, temp_interval, model_motion.cardiac_phase);
        
        % cardiac
        temp_interval_motion(temp_interval_motion > 0) = 0; % distolic
        cMo(time_index_1:time_index_2) = temp_interval_motion;
        
        %     plot(temp_interval_motion)
    end
end

% resp
current_time = 0;
rMo = zeros(size(t_signal));
if model_motion.resp_interval == 0
    % Leave rMo as zeros, no cardiac motion
else
    % Simple sinusoid respiratory motion, with variation
    while (current_time < t_signal(end))
        
        time_index_1 = find(t_signal <= current_time, 1, 'last');
        
        temp_interval = round(model_motion.resp_interval + randn(1)*(model_motion.resp_variation*model_motion.resp_interval/100));
        current_time = current_time + temp_interval;
        
        time_index_2 = find(t_signal <= current_time, 1, 'last');
        t_signal_temp = t_signal(time_index_1:time_index_2) - t_signal(time_index_1);
        
        temp_interval_motion = model_motion.motion_model(t_signal_temp, temp_interval, model_motion.resp_phase);
        
        rMo(time_index_1:time_index_2) = 0.5*temp_interval_motion;
    end
end
model_motion.cMo    = cMo;
model_motion.rMo    = rMo;

%% PLOT Simulation overview
% can move this to sampling function for update
figure,
subplot(2,2,1); hold on;
plot(t_signal, cMo)
plot(t_signal, rMo)
xlabel('Time (ms)'); legend({'cardiac','resp'}); title('Motion model');

subplot(2,2,2); 
if signal_info.ext_sig_flag ==0
    hold on; plot(t_signal, ...
        [tv_signal(tissue_params.fat,       scan_params, t_signal, signal_info.inversion_alpha)' ...
        tv_signal(tissue_params.muscle,    scan_params, t_signal, signal_info.inversion_alpha)' ...
        tv_signal(tissue_params.blood,     scan_params, t_signal, signal_info.inversion_alpha)' ...
        tv_signal(tissue_params.liver,     scan_params, t_signal, signal_info.inversion_alpha)' ...
        tv_signal(tissue_params.scar,      scan_params, t_signal, signal_info.inversion_alpha)']);
else
    plot(t_signal, signal_info.ext_sig_profs);
end

xlabel('Time (ms)'); legend({'fat','muscle','blood','liver','scar'}); title('Signal model');


%% NUFFT set-up
% NUFFT to deal with any trajectory

traj = k;
% == plot trajectory ==
subplot(2,2,3);

hold on, plot(traj(:,:,1), traj(:,:,2), '-','Color', [150, 175, 200]/255), plot(traj(:,:,1), traj(:,:,2), '.','Color', [150, 150, 150]/255)
kr_step = sqrt(2*abs(mean(diff(traj(:,1,1)))));
xlim([min(k(:)) max(k(:))] +kr_step*[-2 2]); ylim([min(k(:)) max(k(:))] +kr_step*[-2 2]); axis image;

plot(traj(:,1,1), zeros( size( traj(:,1,1) )), 'k--')
plot(zeros(size(traj(1,:,2))), traj(1,:,2), 'k--')
plot(0,0,'ko')

plot(traj(:,2,1), traj(:,2,2), 'r-'); % highlight path
title('Trajectory');
% == plot trajectory ==

traj = reshape(traj,samples*shots,2);
traj = traj*(pi/max(traj(:)));

% Prep density compensation weights
ksp = traj./pi.*((enc_info.res/2)./enc_info.FOV);
mask = true(enc_info.res(1),enc_info.res(2));
sizeMask = size(mask);
nufft_args = {sizeMask, [6 6], 2*sizeMask, sizeMask/2, 'table', 2^12, 'minmax:kb'};
G = Gmri(ksp, mask, 'fov', enc_info.FOV, 'basis', {'dirac'}, 'nufft', nufft_args);
wi = abs(mri_density_comp(ksp, 'pipe','G',G.arg.Gnufft));
wi = wi*(sqrt(prod(enc_info.res)))/sum(wi(:));

% Prepare NUFFT operator
st = nufft_init(traj, [enc_info.res], [6 6], [enc_info.res].*2, [enc_info.res]./2);

%% Sample data and add noise

% k needs to be in format of 2 x samples
k = permute( reshape( k, [samples*shots 2]) , [2 1] );

% simple single coil for snr scout
sens.model = 'sinusoidal';
sens.data = 1;
sens.param = 1;

SL = make_SL_cardiac(1, t_signal, enc_info, signal_info, model_motion, tissue_params);

% Sample numerical phantom
m_analytical = MRData(SL,sens,k);
noise_vec = complex(randn(size(m_analytical)),randn(size(m_analytical)));

% Scale image signal (~ pre-whitened)
if isempty(signal_info.snr_scale)
    temp = abs( nufft_adj(m_analytical'.*wi, st) ); % scout transform.. 
    snr_scale = ( ( std(noise_vec(:)) * signal_info.image_snr)/ ( max(temp(:)) ) );
    signal_info.snr_scale = snr_scale;
else
    snr_scale = signal_info.snr_scale;
end

% preview image
    % m_analytical = m_analytical*snr_scale + noise_vec;
    % test_images = abs( nufft_adj(m_analytical'.*wi, st) ); % scout transform.. 
    % figure, imshow(test_images,[0 signal_info.image_snr]);

    %% Multicoil sampling
    subplot(2,2,4);
    
    if enc_info.coils == 0
        % for testing speed, can run without multicoil
        % i.e. using simple single coil deom above
        
        sens_sl{1} = sens;
        NbCoils = 1;
        relative_coil_snr = 1;
        
    else
        % using example, not tested params that much (originally set for brain?)
        sens.model = 'sinusoidal';
        sens.param = 7;
        
        im = RasterizePhantom(SL,enc_info.res,[1],0);
        support = (im>1e-3);numel(find(support));
        NbCoils = enc_info.coils;
        
        if isempty(enc_info.sens_maps)      
            % Generate sensitivity maps
            coil.Nb_coils = enc_info.coils;
            coil.res = enc_info.res;
            coil.type = 'biot';
            coil.param.FS = 0.28;   % FOV width
            coil.param.D = 0.17;    % Distance center->coil
            coil.param.R = 0.05;    % radius of the coil
            %coil.param.rand = 0;   % it doesn't matter what value this is, it will run a rand coil distribution if this is a field!
            coil = simulate_sensitivities(coil);
            % NbCoils = size(coil.sensitivity,3);
            
            sensitivity = coil.sensitivity/max(reshape(abs(coil.sensitivity.*repmat(support,[1,1,NbCoils])),1,numel(coil.sensitivity)));
            enc_info.sens_maps = sensitivity;                   % export
        else
            % Custom sensitivity maps
            sensitivity = enc_info.sens_maps;
        end

        % Inhomogenous coil pickup
        if isempty(enc_info.rel_coil_snr)
            relative_coil_snr = (10 + randn([1 NbCoils]));
            relative_coil_snr = relative_coil_snr./max(relative_coil_snr(:));
            enc_info.rel_coil_snr = relative_coil_snr;          % export
        else % Custom (add correlation?)
            relative_coil_snr = enc_info.rel_coil_snr;
        end
        
        % Generate sensitivities once
        for iCoil = 1:NbCoils
            sens_sl{iCoil} = SensFitting(sensitivity(:,:,iCoil),sens.model,sens.param,support);            
        end
    end
    
%% Run simulation

test_frames = length(t_signal);
test_images = zeros([enc_info.res test_frames]);
kdata = zeros(length(k), test_frames, NbCoils, 'single');

for i = 1:test_frames
    if isunix && ~ismac
    disp(['Current frame ' num2str(i) ' / ' num2str(test_frames)]);
    end
    
    % Call model for signal and motion state
    SL = make_SL_cardiac(i, t_signal, enc_info, signal_info, model_motion, tissue_params) ;
    
    for iCoil = 1:NbCoils
        % Loop through coils for sampling
        if ((i == 1) && (iCoil ==1))
            tic;
        end
        m_analytical_coil = MRData(SL,sens_sl{iCoil},k);
        if ((i == 1) && (iCoil ==1))
            ttt = toc
            disp(['expected time ' num2str(ttt*test_frames*NbCoils/60) ' min /' num2str(ttt*test_frames*NbCoils/3600) ' hrs'])
        end
        m_analytical_coil = relative_coil_snr(iCoil)*m_analytical_coil*snr_scale + complex(randn(size(m_analytical_coil)),randn(size(m_analytical_coil)));
        
        kdata(:,i, iCoil) = m_analytical_coil;
    end
    kdata_temp = squeeze(kdata(:, i, :));
    
    if NbCoils > 1
        test_images(:,:,i) = ismrm_rss(nufft_adj(kdata_temp.*repmat(wi,[1 NbCoils]), st));
    else
        test_images(:,:,i) = abs(nufft_adj(kdata_temp.*repmat(wi,[1 NbCoils]), st));
    end
    imshow(test_images(:,:,i),[0 signal_info.image_snr]); title(['Time ' num2str(t_signal(i)) ' ms']); drawnow;
end

%% package up

kdata = reshape(kdata, [samples shots test_frames NbCoils]);

sim_info.enc_info       = enc_info; 
sim_info.signal_info    = signal_info; 
sim_info.model_motion   = model_motion;
sim_info.test_images    = test_images;

end

function SL = make_SL_cardiac(i, t_signal, enc_info, signal_info, model_motion, tissue_params)
%% Make phantom 

SL.FOV = enc_info.FOV; % the FOV actually changes the shape of the phantom
    %   This head phantom is the same as the Shepp-Logan except the intensities are changed to yield higher contrast in the image.  Taken from Toft, 199-200.
    % fat_signal      = 1; % muscle_signal   = .4; % blood_signal    = .75; % liver_signal    = .3; % scar_signal     = .9;

% unpack
cMo = model_motion.cMo;
rMo = model_motion.rMo;

if signal_info.ext_sig_flag ==0
    %% Apply internal signal evolution
    
    ss_signal       = signal_info.ss_signal;
    tv_signal       = signal_info.tv_signal;
    scan_params     = signal_info.scan_params;
    
    fat_signal      = tv_signal(tissue_params.fat, scan_params, t_signal(i), signal_info.inversion_alpha);
    muscle_signal   = tv_signal(tissue_params.muscle, scan_params, t_signal(i), signal_info.inversion_alpha);
    blood_signal    = tv_signal(tissue_params.blood, scan_params, t_signal(i), signal_info.inversion_alpha);
    liver_signal    = tv_signal(tissue_params.liver, scan_params, t_signal(i), signal_info.inversion_alpha);
    scar_signal     = tv_signal(tissue_params.scar, scan_params, t_signal(i), signal_info.inversion_alpha);
    
else
    %% Extract external signal profiles
    % assuming certain order
    
    fat_signal      = signal_info.ext_sig_profs(i,1);
    muscle_signal   = signal_info.ext_sig_profs(i,2);
    blood_signal    = signal_info.ext_sig_profs(i,3);
    liver_signal    = signal_info.ext_sig_profs(i,4);
    scar_signal     = signal_info.ext_sig_profs(i,5);
    
end

    % overlapping signals
    scar_signal     = scar_signal - muscle_signal; 
    vessel_signal   = blood_signal - liver_signal;
    trab_signal    = muscle_signal - blood_signal;

    %% Apply motion 
    % empirically designed
    %               A            a       b      x0      y0    phi
    %        ---------------------------------------------------------------
shep = [ fat_signal        .90+ rMo(i)*0.0      .75 + rMo(i)*0.070   0                                  0                                   10         % fat layer outer
        -fat_signal        .86+ rMo(i)*0.0      .70 + rMo(i)*0.035   0                                  0                                   15         % fat layer inner
         muscle_signal     .86+ rMo(i)*0.0      .70 + rMo(i)*0.035   0                                  0                                   15         % abdominal muscle layer outer
        -muscle_signal     .79 + rMo(i)*0.0     .63 + rMo(i)*0.07    0                                  0                                   11         % abdominal muscle layer inner
         muscle_signal     .30 + cMo(i)*0.05    .31 + cMo(i)*0.05    0.31 + rMo(i)*0.05                 0.1  + rMo(i)*0.05                  4          % LV outer
        -muscle_signal     .23 + cMo(i)*0.1     .24 + cMo(i)*0.1     0.31 + rMo(i)*0.05                 0.1  + rMo(i)*0.05                  4          % LV inner
         blood_signal      .23 + cMo(i)*0.1     .24 + cMo(i)*0.1     0.31 + rMo(i)*0.05                 0.1  + rMo(i)*0.05                  4          % LV blood
         trab_signal       .05 + cMo(i)*0.03    .07 + cMo(i)*0.05    0.21 - cMo(i)*0.04 + rMo(i)*0.05   0.23 + cMo(i)*0.04 + rMo(i)*0.05   -12         % trabeculae 1
         trab_signal       .052 + cMo(i)*0.028  .068+ cMo(i)*0.034   0.37 + cMo(i)*0.03 + rMo(i)*0.05  -0.02 - cMo(i)*0.02 + rMo(i)*0.05    23         % trabeculae 2
         scar_signal       .03                  .03                  0.47 + cMo(i)*0.04 + rMo(i)*0.05   0.32 + cMo(i)*0.035 + rMo(i)*0.05   4          % scar
         liver_signal      .34 + rMo(i)*0.03    .45 + rMo(i)*0.04   -0.34 - rMo(i)*0.05                -0.1  + rMo(i)*0.05                  33         % liver
         vessel_signal     .03                  .14                 -0.31 + rMo(i)*0.05                -0.2  - rMo(i)*0.05                  62         % Liver vessel 1
         vessel_signal     .01                  .08                 -0.27 + rMo(i)*0.05                -0.1  - rMo(i)*0.05                 -21         % Liver vessel 2
        ];
    
    %% build phantom structure
    SL.region = cell(1,size(shep,1));
    for ii = 1:size(shep,1)
        SL.region{ii} = struct('type','ellipse',...
            'center',[-shep(ii,5),shep(ii,4)]/2,... % in FOV units
            'angle',shep(ii,end)*pi/180,...
            'weight',shep(ii,1),...
            'width',[shep(ii,[3,2])]); % in FOV units
    end
    clear shep ii;

end

function input_help

description_cell = {
'\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
' '
'<strong>sim_SL_cardiac</strong>'
' '
'The cardiac model was a simple attempt at just using ellipses '
'Motion is modelled by sinusoid behaviour for both cardiac and respiratory'
'Cardiac motion is a "one-directional" sine (just -ve, 0 for +ve)'
'The tissue parameters are coded in to the simulator'
' '
'function [<strong>kdata</strong> [samples x shots x numAcqs x coils], <strong>sim_info</strong> [struct]] ...'
' = sim_SL_cardiac(<strong>enc_info</strong>, <strong>signal_info</strong>, <strong>model_motion</strong>)'
' '
'\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
'Choose from: '
' '
'<strong>sampling design</strong>' 
'enc_info.res            = [256 256];'
'enc_info.FOV            = [32 32];'
'enc_info.k              = GenerateFullCart2DKspace(enc_info.res,enc_info.FOV); % Arbitrary trajectory (everything nuffted)'
'enc_info.k              = reshape(permute(enc_info.k,[2 1]),[enc_info.res 2]); % nufft format : [samples shots 2]'
'enc_info.coils          = 6;'
'enc_info.rel_coil_snr   = [];                          % defined relative coil SNR'
'enc_info.sens_maps      = [];'
' '
'<strong>contrast</strong>' 
'signal_info.TR              = 10;'
'signal_info.numAcqs         = 24;'
'signal_info.contrast        = 1;'
'signal_info.time_var_sig    = 0;'
'signal_info.inversion_alpha = 1.8;'
'signal_info.image_snr       = 50;'
'signal_info.snr_scale       = [];                      % custom/consistent scale factor to transform'                     
'signal_info.ext_sig_flag    = 0;                       % [0]: simple model, [1]: arbitrary signal model'                      
'signal_info.ext_sig_profs   = [];                      % [NumAcqs x 5] : [fat muscle blood liver scar]'
' '
'<strong>motion</strong>' 
'model_motion.cardiac_interval   = 777;                 % ms'
'model_motion.cardiac_variation  = 6;                   % percentage'
'model_motion.cardiac_phase      = randn(1);            % offset to starting phase'
'model_motion.resp_interval      = 3318;                % ms'
'model_motion.resp_variation     = 13;                  % percentage'
'model_motion.resp_phase         = randn(1);            % offset to starting phase'
'model_motion.motion_model       = @(t_ms, int_ms, phase) sin(2*pi*((t_ms/int_ms) + phase));'
' ' 
['Click <a href="matlab:matlab.desktop.editor.openAndGoToLine(''' mfilename('fullpath') '.m'',560);">here</a> for dependencies']
' '
'To run defaults:' 
'enc_info.dummy      = 1;'
'signal_info.dummy   = 1;'
'model_motion.dummy  = 1;'
'[kdata, sim_info]  = sim_SL_cardiac(enc_info, signal_info, model_motion)'
};

fprintf(1, '%s \n', description_cell{:});
fprintf(1, '\n');
disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');

end

function dep_help
disp(1);

%%REQUIREMENTS
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
% add all paths:
    % mriphantom_path = uigetdir(matlabroot, 'Select ISMRM sunrise dir');
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

end