function [x] = dicom_read_shim(path1)
% function [x, xinfo] = dicom_read_shim(<path1>)
% Grab the shim settings from the dicom.
% Written for Siemens Aera VE11C
%
% [path1] <optional> to sorted folder containing dicom image
% Otherwise, UIGETDIR will pop-up
%
% First Order shims location:
% [sGRADSPEC.asGPAData0.lOffsetX sGRADSPEC.asGPAData0.lOffsetY sGRADSPEC.asGPAData0.lOffsetZ]
% Second Order shims location:
% [sGRADSPEC.alShimCurrent[0]...sGRADSPEC.alShimCurrent[4]]
%
% Notes
% [1] To sort a folder of dicoms in to scans: > dicom_sort_folder
% [2] Shim values are reported 10x that in the Adjustments UI (units?)
% [3] It has been observed that an alShimCurrent value may be skipped..
%
% R Ramasawmy Oct 2019

%% Choose data
if nargin < 1
    path1 = uigetdir;
end

% Load first image in series // can switch this to picking specfic file?
% Just makes it easier when writing a script..
dir1 = dir(path1); temp = 0; i = 0;
while temp < 1
    i = i + 1;
    if ~isempty(regexpi(dir1(i).name, '.IMA')) || ~isempty(regexpi(dir1(i).name, '.DCM'))
        %         x = dicomread([path1 filesep dir1(i).name]);
        %         figure, imshow(x,[]);
        disp(' '); disp('================================'); disp(' ');
        disp(['Reading ' path1 ' -image 1']);
        fileID = fopen([path1 filesep dir1(i).name]);
        dicomraw = fscanf(fileID, '%s');
        
        temp = 1;
    end
end
fclose(fileID); clear fileID temp dir path1; 

%% First Order Shims

a = regexp(dicomraw, 'sGRADSPEC.asGPAData');

for i = 1:length(a)
    temp = dicomraw(a(i):(a(i) + 90)); % heuristic
    
    if ~isempty(regexp(temp, 'lOffset')) % remove to grab everything
    
    b = regexp(temp, '=');
    if ~isempty(b)
        
        numpos = b(1)+1;
        
        try
            next_one = regexp(temp(numpos:end), 's'); % or sGRADSPEC sTXSPEC
            next_one = next_one(1);
            svalue_x = temp(numpos:(next_one + numpos - 2));
            value_x = str2double(temp(numpos:(next_one + numpos - 2)));
        catch
            c = find(isstrprop(temp(numpos:end), 'digit')==0, 1, 'first');
            svalue_x = temp(numpos:(numpos+c-2));
            value_x = str2double(temp(numpos:(numpos+c-2)));
        end
        
        name_x = temp(1:b-1);
        d = regexp(name_x, '_');
        name_x(d)  = [];
        d = regexp(name_x, '[');
        name_x(d)  = [];
        d = regexp(name_x, ']');
        name_x(d)  = [];
        
        
        eval([name_x '=' svalue_x '*0.1;']);
    end
    end
end

% Not yet observed in linear shims, but add a catch for missing values.
shim_dir = ['X', 'Y', 'Z'];
for i = 1:3
    if ~isfield(sGRADSPEC.asGPAData0, ['lOffset' shim_dir(i)])
        sGRADSPEC.asGPAData0.(matlab.lang.makeValidName(['lOffset' shim_dir(i)])) = 0;
        warning(['asGPAData0.lOffset' shim_dir(i) ' does not have a value! Setting to 0']);
    end
end

% [sGRADSPEC.asGPAData0.lOffsetX sGRADSPEC.asGPAData0.lOffsetY sGRADSPEC.asGPAData0.lOffsetZ]

%% Second Order Shims

a = regexp(dicomraw, 'sGRADSPEC.alShimCurrent');

for i = 1:length(a)
    temp = dicomraw(a(i):(a(i) + 60)); % heuristic
    b = regexp(temp, '=');
      if ~isempty(b)
    numpos = b(1)+1;
    
    try
        next_one = regexp(temp(numpos:end), 's'); % or sGRADSPEC sTXSPEC
        next_one = next_one(1);
        svalue_x = temp(numpos:(next_one + numpos - 2));
        value_x = str2double(temp(numpos:(next_one + numpos - 2)));
    catch
        c = find(isstrprop(temp(numpos:end), 'digit')==0, 1, 'first');
        svalue_x = temp(numpos:(numpos+c-2));
        value_x = str2double(temp(numpos:(numpos+c-2)));
    end
    
    name_x = temp(1:b-1);
    d = regexp(name_x, '_');
    name_x(d)  = [];
    d = regexp(name_x, '[');
    name_x(d)  = [];
    d = regexp(name_x, ']');
    name_x(d)  = [];
    
    
    eval([name_x '=' svalue_x '*0.1;']);
      end
end

% fill in missing values? (Why does this happen!?)
for i = 0:4
    if ~isfield(sGRADSPEC, ['alShimCurrent' num2str(i)])
        sGRADSPEC.(matlab.lang.makeValidName(['alShimCurrent' num2str(i)])) = 0;
        warning(['alShimCurrent' num2str(i) ' does not have a value! Setting to 0']);
    end
end

% [sGRADSPEC.alShimCurrent0 sGRADSPEC.alShimCurrent1 sGRADSPEC.alShimCurrent2 sGRADSPEC.alShimCurrent3 sGRADSPEC.alShimCurrent4]

%% Print to screen
% clc
disp(' '); disp('================================'); disp(' ');
Shim_Labels = {'lOffsetX', 'lOffsetY', 'lOffsetZ', 'alShimCurrent0', 'alShimCurrent1', 'alShimCurrent2', 'alShimCurrent3', 'alShimCurrent4'}';
Value = [sGRADSPEC.asGPAData0.lOffsetX sGRADSPEC.asGPAData0.lOffsetY sGRADSPEC.asGPAData0.lOffsetZ sGRADSPEC.alShimCurrent0 sGRADSPEC.alShimCurrent1 sGRADSPEC.alShimCurrent2 sGRADSPEC.alShimCurrent3 sGRADSPEC.alShimCurrent4]';
disp(table( Shim_Labels,Value ));
disp(' '); disp('================================'); disp(' ');

% if exist('printstruct', 'file')
%     printstruct(sGRADSPEC)
% end

%% Return values
x = Value;

end