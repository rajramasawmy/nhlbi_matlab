function [x, xinfo] = dicom_load_scan(path1, scanRange)
% [x, xinfo] = load_scan(<optional-UI>path1, <optional>scanRange)
%
% Example Usage: 
%     x = dicom_load_scan; > use UI
%     [x, info] = dicom_load_scan('directory of dicom scan'); 
%
% Returns a 3D array <x>, and the dicom header <xinfo>, given the directory
% of a folder containing DICOM images. (Requires all the images to have the
% same dimensions.)
% Raj Ramasawmy, NHBLI 2017

if nargin < 1
    path1 = uigetdir;
    dir1 = dir(path1);
    NumObj = length(dir1) - 2;
    scanRange = 1:NumObj;
elseif nargin < 2
    dir1 = dir(path1);
    NumObj = length(dir1) - 2;
    scanRange = 1:NumObj;
end
dir1 = dir(path1);

info_flag = 0; te_vec = 0;
for j = scanRange
    % Different layers can be extracted:
    %     [x(:,:,j), a(:,:,j),  b(:,:,j),  c(:,:,j)] = dicomread([path1 filesep dir1(j+2).name]);
    x(:,:,j) = dicomread([path1 filesep dir1(j+2).name]);

    if ~info_flag
        % This bit is 1) very specific for MRI, and for
        % multi-echo/-inversion time experiments. 
        xinfo = dicominfo([path1 filesep dir1(j+2).name]);
        
        if isfield(xinfo, 'EchoTime'); if (~isempty(xinfo.EchoTime)); te_vec(j) = xinfo.EchoTime; end; end

        if isfield(xinfo, 'InversionTime')
            
            % extract inversion times >> this could be generalised to a
            % string input as I rarely use the "Scan Range" function
            temp(j) = xinfo.InversionTime;
        else
            info_flag = 1;
        end
        
    end
end

if exist('xinfo', 'var')
    if isfield(xinfo, 'InversionTime')
    xinfo.InversionTime = temp;
    end
end
if exist('te_vec', 'var')
if(length(unique(te_vec))==1)
    xinfo.EchoTime = te_vec;
end
end

% === patient info === 
% \\\ Anonymize header
xinfo.PatientName = '';
xinfo.PatientBirthDate = '';
xinfo.PatientID = '';
xinfo.Filename = '';
xinfo.InstitutionAddress = '';
xinfo.ReferringPhysicianName = '';
% \\\ Keeping this data;
% xinfo.PatientSex
% xinfo.PatientAge
% xinfo.PatientSize
% xinfo.PatientWeight
% xinfo.PatientPosition

% === increase precision for calculations ===
x = double(x);

end

