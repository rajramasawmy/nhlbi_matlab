function [img_s] = recon_spiral(dfile,  nfile)
% function [data_struct] = recon_spiral(data_file,  noise_file)
% 
% Simple spiral recon using MRD data format 
%
% Requirements: 
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
% 
% ------------------------------
% The ISMRMRD matlab class
% ------------------------------
    % from the MRD/ISMRMRD community
%
% available at
% https://github.com/ismrmrd/ismrmrd/tree/master/matlab
% And to add:
    % ismrmrd_path = uigetdir(matlabroot, 'Select ismrmrd-matlab dir');
    %
    %     addpath([ismrmrd_path])
% 
%
% R Ramasawmy NHLBI June 2020
% 

%% Set up

make_nhlbi_utils;

% ============================================
% Load data
% ============================================

dsetin = ismrmrd.Dataset(dfile, 'dataset');
MRD_h = ismrmrd.xml.deserialize(dsetin.readxml);
disp(['Reconstructing: ' MRD_h.measurementInformation.protocolName]);

raw_data    = dsetin.readAcquisition; % h5read(dfile,'/dataset/data');

% nhlbi_utils.plot_experiment(raw_data);

% ============================================
% Grab imaging parameters
% ============================================

interleaves     = 1 + (MRD_h.encoding.encodingLimits.kspace_encoding_step_1.maximum);
pe2             = 1 + (MRD_h.encoding.encodingLimits.kspace_encoding_step_2.maximum);
averages        = 1 + (MRD_h.encoding.encodingLimits.average.maximum);
slices          = 1 + (MRD_h.encoding.encodingLimits.slice.maximum);
contrasts       = 1 + (MRD_h.encoding.encodingLimits.contrast.maximum);
phases          = 1 + (MRD_h.encoding.encodingLimits.phase.maximum);
sets            = 1 + (MRD_h.encoding.encodingLimits.set.maximum);
reps            = 1 + (MRD_h.encoding.encodingLimits.repetition.maximum);

samples         =      double(raw_data.head.number_of_samples(1));
dt              =      raw_data.head.sample_time_us(1)*1e-6;
channels        =      double(raw_data.head.active_channels(1));

matrix          = MRD_h.encoding.reconSpace.matrixSize.x;                   % assuming a square spiral matrix
matrix_size     = [matrix matrix];

disp(' ');disp('### Experiment Dimensions ###');disp(' ');
Experiment_parameters = {'Samples', 'Interleaves', 'PE2', 'Averages', 'Slices', 'Contrasts', 'Phases', 'Repetitions', 'Sets', 'Channels'}';
Value = [samples interleaves pe2 averages slices contrasts phases reps sets channels]';
disp(table( Experiment_parameters,Value )); clear Experiment_parameters Value; disp(' ');

%% Noise checks

if isempty(nfile)
    dmtx = diag(ones(1,channels));
else
    dmtx = nhlbi_utils.noise_adjust(nfile, MRD_h, dt);
end

%% Build Nominal Fully Sampled traj and gradients

FOV = MRD_h.encoding.reconSpace.fieldOfView_mm.x/10;

traj_setup.gMax = MRD_h.encoding.trajectoryDescription.userParameterDouble(1).value;
traj_setup.sMax = MRD_h.encoding.trajectoryDescription.userParameterDouble(2).value;

krmax = 1/(2*(FOV(1)/matrix_size(1)));

[k,g] = vds(traj_setup.sMax, traj_setup.gMax, dt, interleaves, FOV, krmax); close;


%% Rotate spiral design 
% crop to data
if samples > length(k)
    samples2 = length(k);
else
    samples2 = samples;
end
trajectory_nominal = zeros(samples2,interleaves,2);
gradients_nominal =  zeros(samples2,interleaves,2);

neg = 1;
for solid_int= 1:interleaves
    rot = (solid_int-1)*(2*pi/interleaves);
    trajectory_nominal(:,solid_int,2) = neg*-( real(k(1:samples2)) *cos(rot) + imag(k(1:samples2)) *sin(rot));
    trajectory_nominal(:,solid_int,1) = neg*-(-real(k(1:samples2)) *sin(rot) + imag(k(1:samples2)) *cos(rot));
    gradients_nominal(:,solid_int,2)  = neg*-( real(g(1:samples2)) *cos(rot) + imag(g(1:samples2)) *sin(rot));
    gradients_nominal(:,solid_int,1)  = neg*-(-real(g(1:samples2)) *sin(rot) + imag(g(1:samples2)) *cos(rot));
end

%% GIRF corrections

if ~(exist('apply_GIRF', 'file')==0)
   
    % ============================================
    % Apply GIRF
    % ============================================
    
    R_pcs_gcs = [raw_data.head.phase_dir(:,1), raw_data.head.read_dir(:,1), raw_data.head.slice_dir(:,1)  ]; %Rotation matrix
   
    sR.R = R_pcs_gcs;
    sR.T = MRD_h.acquisitionSystemInformation.systemFieldStrength_T;
    
    trajectory_nominal = apply_GIRF(gradients_nominal, dt, sR); 
    trajectory_nominal = trajectory_nominal(:,:,1:2);
    
end

%% Collect all data
kspace = complex(zeros([samples interleaves pe2 averages slices contrasts phases reps sets channels],'single'));
% disp(['Kspace dims: ' num2str(size(kspace))])

for ii = 1:length(raw_data.data)
    
    d1 = raw_data.data{ii};
    d1 = ismrm_apply_noise_decorrelation_mtx(d1, dmtx);
    
    kspace(:,...
        raw_data.head.idx.kspace_encode_step_1(ii)+1, ...
        raw_data.head.idx.kspace_encode_step_2(ii)+1, ...
        raw_data.head.idx.average(ii)+1, ....
        raw_data.head.idx.slice(ii)+1 , ...
        raw_data.head.idx.contrast(ii)+1, ...
        raw_data.head.idx.phase(ii)+1, ...
        raw_data.head.idx.repetition(ii)+1, ...
        raw_data.head.idx.set(ii)+1, ...
        :) = d1;
    
end

%% recon

% ============================================
% NUFFT operator set-up
% ============================================

% Assuming same encoding per image!
omega = trajectory_nominal*(pi/max(max(max(trajectory_nominal))));
omega = reshape(omega,interleaves*samples2,2);
recon_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);

% ============================================
% spiral-out weight estimation
% ============================================
traj_info.traj = trajectory_nominal;
traj_info.grad = gradients_nominal;
recon_weights = nhlbi_utils.ismrm_calculate_weights(traj_info, 'Hoge');

kspace = mean(kspace,4); % pseudo-rep will need to preceed this step

% ============================================
% Main reconstruction
% ============================================

figure,
cCoil_imgs = zeros([matrix_size pe2 1 slices contrasts phases reps sets ]);

if pe2 > 1 
    % stack of spirals recon
    kspace = fftshift( ifft( ifftshift(kspace), pe2, 3) );
end
    
for par = 1:pe2
    for slc = 1:slices
        for coc = 1:contrasts
            for phc = 1:phases
                for repc = 1:reps
                    for setc = 1:sets
                        
                        % ============================================
                        % Gridding and B1-coil combination
                        % ============================================
                        
                        data_temp = squeeze(kspace(1:samples2,:,par,1,slc,coc,phc,repc,setc,:));
                        
                        data_temp = reshape(data_temp,interleaves*samples2,channels);
                        x = nufft_adj(data_temp.*repmat(recon_weights,[1,channels]), recon_st)*sqrt(averages);
                        img_coil = squeeze(x);
                        
                        csm = ismrm_estimate_csm_walsh( img_coil );
                        ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened
                        cCoil_imgs(:,:,par,1,slc,coc,phc,repc,setc)= abs( sum( squeeze( img_coil ) .* ccm_roemer_optimal, 3) );

                        % ============================================
                        % View updated figure
                        % ============================================

                        imshow(quick_scale(cCoil_imgs(:,:,par,1,slc,coc,phc,repc,setc)),[0 4]); 
                        drawstring = [];
                        if pe2 > 1
                            drawstring = [drawstring 'Slice ' num2str(par) ' '];
                        end
                        if slices > 1
                            drawstring = [drawstring 'Slice ' num2str(slc) ' '];
                        end
                        if contrasts > 1
                            drawstring = [drawstring 'Contrast ' num2str(coc) ' '];
                        end
                        if phases > 1
                            drawstring = [drawstring 'Phase ' num2str(phc) ' '];
                        end
                        if reps > 1
                            drawstring = [drawstring 'Repetition ' num2str(repc) ' '];
                        end
                        if sets > 1
                            drawstring = [drawstring 'Set ' num2str(setc) ' '];
                        end
                        title(drawstring); drawnow;
                                                
                    end
                end
            end
        end
    end
end

imgs = squeeze(cCoil_imgs); close; 

%% return components
img_s.img = imgs;
img_s.header = MRD_h;


end

function y = quick_scale(x)
y = x./mean(x(:));
end