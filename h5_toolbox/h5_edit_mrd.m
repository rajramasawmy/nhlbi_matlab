function [sfile] = h5_edit_mrd(dfile, user_opts)
% [sfile] = h5_edit_mrd(<dfile, UI>, <user_opts, UI>)
% 
% *** THIS FILE IS JUST SOME VERY SPECIFIC EXAMPLES ***
%
% 20200429 : Current implementation only automatically loads trajectories
%
% ///////////////////////////////////////////////////////////////////////
% This relies on the user prepping the trajectories accordingly!!
%
% whos ExampleTrajFile.mat
%       <strong>traj_data</strong> : [numAcquisitions samples*traj_dims] 
%           The trajectory data must be:
%                 - ordered for the scan 
%                   i.e. a dataset with 2-repetitions, 8-shots, with 2000 samples and w/ weights will have dims:
%                   [8*2 2000*3] 
%                 - |max val|=0.5
%                 - sorted such that the x/y/weights are in the inner 
%                   loop of the sample, 
%                   i.e. the 1st interleave x-traj = traj_data(1:traj_dims:end)
%       <strong>traj_dims</strong>: 2 (2D, no weights) / 3 (2D with weights)
%
% If the spiral is designed with vds, can also use header info automatically
%
% ///////////////////////////////////////////////////////////////////////
% ??? Might be easiest to just stop program and let the user do whatever?
% ///////////////////////////////////////////////////////////////////////
%
% example use:
% >> h5_edit_mrd                % Use UI for file load (load file to edit, load traj file)
%
% % defaults user_opts
% % user_opts.header        = [];
% % user_opts.data          = [];
% % user_opts.traj          = 1;
% % user_opts.traj_and_sort = 0; % being too lazy to do it yourself <2, 3>
% % user_opts.waveform      = 0;
% % user_opts.savename      = []; % 'out_edit.h5';
%
% VDS example (i.e. offline add GIRF)
% >> user_opts.traj = 'vds';
% Add output filename
% >> user_opts.savename = 'OUTPUT_FILE.h5';
% h5_edit_mrd('INPUT_FILE.h5', user_opts);
% 
% using ismrmrd/matlab methods
% https://github.com/ismrmrd/ismrmrd > https://github.com/ismrmrd/ismrmrd/tree/master/matlab/%2Bismrmrd
%
% ///////////////////////////////////////////////////////////////////////
% 
% R Ramasawmy, NHLBI

%% Add ismrmrd-matlab methods
if exist('ismrmrd.Dataset')~=8
    rLoc = mfilename('fullpath');
    temp = regexp(rLoc, filesep);
    rLoc = [rLoc(1:temp(end))];
    addpath(genpath(rLoc))
end


%% Load file to edit

if nargin < 1
    [fname, dirPath] = uigetfile('*.*', 'Choose data .h5 file');
    dfile = [dirPath fname];
    
% %     % dialogue box version
% %     list = {'Attach traj','Attach waveform','Edit header'};
% %     [indx,tf] = listdlg('ListString',list);
% %     if tf == 0
% %         user_opts_in.traj = 0;
% %         user_opts_in.waveform = 0;
% %         user_opts_in.header = 0;
% %     else
% %         user_opts_in.traj = sum(ismember(indx,1));
% %         user_opts_in.waveform = sum(ismember(indx,2));
% %         user_opts_in.header = sum(ismember(indx,3));
% %     end
end

if exist('user_opts','var')
    user_opts_in = user_opts;
end    
user_opts_in_list = fieldnames(user_opts_in);

% defaults
user_opts.header        = [];
user_opts.data          = [];
user_opts.traj          = 1;
user_opts.traj_and_sort = 0; % being too lazy to do it yourself <2, 3>
user_opts.waveform      = 0;
user_opts.savename      = []; % 'out_edit.h5';

% overwrite with options 
for i = 1:length(user_opts_in_list)
    if isfield(user_opts, user_opts_in_list{i})
        user_opts.(matlab.lang.makeValidName(user_opts_in_list{i})) = user_opts_in.(matlab.lang.makeValidName(user_opts_in_list{i}));
    end
end

clear user_opts_in*

disp('<strong>User options</strong>');
disp(user_opts); disp(' ');

%% Load data

disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');
disp('Input h5');
h5disp(dfile);
disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');

dfile_info  = h5info(dfile);
dfile_dsets = dfile_info.Groups.Datasets;

dsetin = ismrmrd.Dataset(dfile, 'dataset');
hdr = ismrmrd.xml.deserialize(dsetin.readxml);
xmlstring = ismrmrd.xml.serialize(hdr); % xmlstring   = h5read(dfile,'/dataset/xml');

raw_data    = h5read(dfile,'/dataset/data'); % D = dsetin.readAcquisition();

if length(dfile_dsets) == 3 % some assumptions here
wav_data    = h5read(dfile,'/dataset/waveforms');
end

%% edit header

if (~isempty(user_opts.header))
    disp('EDITING HEADER'); 
    
    if isstruct((user_opts.header))
        % If only editing specific header elements
    
        hdr2 = user_opts.header;
        list1 = fieldnames(user_opts.header);
       
        for i = 1:length(list1)
            subhdr = update_header(hdr.(matlab.lang.makeValidName(list1{i})), hdr2.(matlab.lang.makeValidName(list1{i})));
            hdr.(matlab.lang.makeValidName(list1{i})) = subhdr;
        end
        
    elseif (ischar(user_opts.header) && regexp(user_opts.header, '.mat'))
        % Load pre-edited header (all of it!)
        
        loadname = whos('-file',user_opts.header);
        eval(['hdr = ' loadname.name]);
    else
        % edit header as you like ??
        
        openvar('hdr');pause; % ???
    end

    xmlstring = ismrmrd.xml.serialize(hdr);
end

%% edit data
if (~isempty(user_opts.data))
disp('EDITING DATA'); 
% some hard coded instructions..
% i.e. crop user_opts.data = [low high]
new_samples = length(user_opts.data(1):user_opts.data(2));
channels = single(raw_data.head.active_channels(1));
samples = single(raw_data.head.number_of_samples(1));
for i = 1:length(raw_data.data)
    d = complex(raw_data.data{i}(1:2:end),raw_data.data{i}(2:2:end));
    d = reshape(d,[samples, channels]);
    d_crop = d(user_opts.data(1):user_opts.data(2),:);
    d_crop = reshape(d_crop, [1 new_samples*channels]);
    d_crop = cat(1, real(d_crop), imag(d_crop));
    d_crop = reshape(d_crop, [2*new_samples*channels 1]);
    raw_data.data{i} = d_crop;  
    raw_data.head.number_of_samples(i) = uint16(new_samples);
end

end

%% edit traj
% ============================================
% TRAJ
% ============================================
    % Note, trajectory is not generic
    % ISMRMRD test "simple_spiral.h5" has dimensions [2000x1 single] for
    % 1000 samples {kx = 1:2:end, ky = 2:2:end} - no weights are attached.
    % Modern Gadgetron (Apr 2020) doesn't have 3D spiral implemented (yet),
    % so the 3rd dimension is assigned to the weights.

if ischar(user_opts.traj) % user_opts.traj == 'vds'
    disp('EDITING TRAJ'); 
    if (~isempty(regexp(user_opts.traj, '.mat')))
        load(user_opts.traj);
        
        if ~exist('traj_dims')
            traj_dims = size(traj_data,2)./raw_data.head.number_of_samples(1);
        end
    
    else
        % override with demo
        % No weights (2D traj dim = 2)
        % [traj_data, traj_dims] = demo_traj(xmlstring, raw_data);
        % add weights (2D traj dim = 3)
        [traj_data, traj_dims] = demo_traj(hdr, raw_data, true);
    end
else
    disp('EDITING TRAJ'); 
    if numel(user_opts.traj) == 1
    % user input
    [traj_file, temp] = uigetfile('*.mat', 'Choose traj mat file, cancel for vds');
    
    if ischar(traj_file)
        load([temp traj_file]); clear traj_file temp; 
    else
        [traj_data, traj_dims] = demo_traj(hdr, raw_data);
    end
    
    if ~exist('traj_dims')
        traj_dims = size(traj_data,2)./raw_data.head.number_of_samples(1);
    end
    
    else
        % Lazily attaching the trajectories here
        traj_dims = user_opts.traj_and_sort;
        traj_data = user_opts.traj*(0.5/max(max(max(user_opts.traj)))); % GT wants -0.5/0.5
        
        [samples, interleaves,~] = size(traj_data);
        
        matrix          = hdr.encoding.reconSpace.matrixSize.x;
        matrix_size     = [matrix matrix];
        FOV             = hdr.encoding.reconSpace.fieldOfView_mm.x/10;
        
        if user_opts.traj_and_sort == 3
            [temp,traj_dims] = add_weights(traj_data, interleaves, samples, matrix_size, FOV);
        end
        traj_data = format_traj_data(temp, traj_dims, samples, interleaves, raw_data);
        
    end
end

for i = 1:length(raw_data.data)
    % traj_data has been arranged so the interleave index matches acq order
    raw_data.traj{i} = squeeze(traj_data(i,:))';
    
    % trajectory "dimensions": https://github.com/gadgetron/gadgetron/blob/master/gadgets/spiral/SpiralToGenericGadget.cpp
    % traj_dims = 2 (2D, no weights), 3 (2D, w/ weights)
    raw_data.head.trajectory_dimensions(i) = uint16(traj_dims); 
end

%% edit waveform

% NO EXAMPLE YET

%% save data

disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');
disp('<strong>Saving h5</strong>');
if isempty(user_opts.savename)
    [sfile, spath] = uiputfile('*.h5');
    sfile = [spath sfile];
else
    sfile = user_opts.savename;
end
disp(['to: ' sfile]); disp(' ');

% Initiate dataset 
dset = ismrmrd.Dataset(sfile);

% xml
disp('Writing <strong>xml</strong>');
dset.writexml(xmlstring);

% data
disp('Appending <strong>acquisition</strong>');     % temp_acq =ismrmrd.Acquisition(raw_data.head, raw_data.traj, raw_data.data);
dset.appendAcquisition(ismrmrd.Acquisition(raw_data.head, raw_data.traj, raw_data.data));

% waveform
if length(dfile_dsets) == 3 
disp('Appending <strong>waveform</strong>');
% untested:
dset.appendWaveform(ismrmrd.Waveform(wav_data.head, wav_data.data));
end

% Write the dataset
disp('Writing h5 file!');
dset.close();

% check
% h5disp(sfile)


end

%% local functions

function [subhdr3] = update_header(subhdr, subhdr2)
% ============================================
% Update specific fields
% ============================================   

subhdr3 = subhdr;
if isstruct(subhdr2) && (length(subhdr3) < 2)
    % recursive search for exising fields
    list1 = fieldnames(subhdr2);
    for i = 1:length(list1) 
        subhdr3.(matlab.lang.makeValidName(list1{i})) = ...
            update_header(subhdr3.(matlab.lang.makeValidName(list1{i})), ...
                          subhdr2.(matlab.lang.makeValidName(list1{i})));
    end
else
    % update parameter
    subhdr3 = subhdr2;
end
end

function [traj_data, traj_dims] = demo_traj(hdr, raw_data, addWeights)
% ============================================
% Get params
% ============================================   

    if nargin < 3
        addWeights = false;
    end
    
    interleaves     = (1 + double(max(raw_data.head.idx.kspace_encode_step_1)));
    samples         =      double(raw_data.head.number_of_samples(1));
    matrix          = hdr.encoding.reconSpace.matrixSize.x;
    dt              = raw_data.head.sample_time_us(1)*1e-6;
    matrix_size     = [matrix matrix];

    FOV = hdr.encoding.reconSpace.fieldOfView_mm.x/10;
    traj_setup.gMax = hdr.encoding.trajectoryDescription.userParameterDouble(1).value;
    traj_setup.sMax = hdr.encoding.trajectoryDescription.userParameterDouble(2).value;%hdr.encoding.trajectoryDescription.userParameterLong.value;
    
    krmax = 1/(2*(FOV(1)/matrix_size(1)));
    
    traj_dims = 2;

% ============================================
% Make traj
% ============================================

    [k,g] = vds(traj_setup.sMax, traj_setup.gMax, dt, interleaves, FOV, krmax); close;
    
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
        trajectory_nominal(:,solid_int,1) = neg*-( real(k(1:samples2)) *cos(rot) + imag(k(1:samples2)) *sin(rot));
        trajectory_nominal(:,solid_int,2) = neg*-(-real(k(1:samples2)) *sin(rot) + imag(k(1:samples2)) *cos(rot));
        gradients_nominal(:,solid_int,1)  = neg*-( real(g(1:samples2)) *cos(rot) + imag(g(1:samples2)) *sin(rot));
        gradients_nominal(:,solid_int,2)  = neg*-(-real(g(1:samples2)) *sin(rot) + imag(g(1:samples2)) *cos(rot));
    end

% ============================================
% Apply GIRF
% ============================================

tRR = 0; % custom clock-shift
R = [raw_data.head.read_dir(:,1),raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ]; %Rotation matrix     
sR.R = R;
sR.T = hdr.acquisitionSystemInformation.systemFieldStrength_T;

trajectory_nominal = apply_GIRF(gradients_nominal, dt, sR, tRR );           % legacy > % trajectory_nominal = apply_GIRF(gradients_nominal, dt, R, tRR );
trajectory_nominal = trajectory_nominal(:,:,1:2);
trajectory_nominal = trajectory_nominal*(0.5/max(max(max(trajectory_nominal)))); % GT wants -0.5/0.5

temp = trajectory_nominal;

if addWeights
    [temp,traj_dims] = add_weights(trajectory_nominal, interleaves, samples2, matrix_size, FOV);
end

traj_data = format_traj_data(temp, traj_dims, samples2, interleaves, raw_data);

end

function [temp,traj_dims] = add_weights(trajectory_nominal, interleaves, samples, matrix_size, FOV)
 
    % ============================================
    % Calc weights
    % ============================================
    
    traj_dims = 3;
    
    omega = trajectory_nominal*2;                                        % local gridder wants kr(1/cm)
    omega = reshape(omega,interleaves*samples,2);
    
    ksp = omega.*((matrix_size/2)./[FOV(1) FOV(1)]);
    mask = true(matrix_size);
    sizeMask = size(mask);
    nufft_args = {sizeMask, [6 6], 2*sizeMask, sizeMask/2, 'table', 2^12, 'minmax:kb'};
    G = Gmri(ksp, mask, 'fov', FOV, 'basis', {'dirac'}, 'nufft', nufft_args);
    
    recon_weights = abs(mri_density_comp(ksp, 'pipe','G',G.arg.Gnufft));
    
    % scale weights such that images ~= SNR-scaled
    recon_weights = recon_weights*sqrt(prod(matrix_size))/(sum(recon_weights(:)));
    
    recon_weights = reshape(recon_weights, samples, interleaves);
    % figure, plot(recon_weights)
    
    temp = cat(3,trajectory_nominal, recon_weights);
    
end

function traj_data = format_traj_data(temp, traj_dims, samples, interleaves, raw_data)
% ============================================
% Send it back
% ============================================

temp = permute(temp,[3 1 2]);
temp = reshape(temp, [traj_dims*samples, interleaves]);
temp = single(temp); % h5 format

%      figure, subplot(2,2,1); plot(temp(:,1))
%      for i = 1:traj_dims
%      subplot(2,2,i+1); hold on, plot(temp(i:traj_dims:end,1));plot(temp(i:traj_dims:end,floor(interleaves/2)));
%      end

traj_data = zeros(length(raw_data.data),size(temp,1), 'single');
for i = 1:length(raw_data.data)
    temp_interleave = 1+raw_data.head.idx.kspace_encode_step_1(i);
    traj_data(i, :) = temp(:,temp_interleave);
end

end