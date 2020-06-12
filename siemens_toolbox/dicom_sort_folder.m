function dicom_sort_folder(path1)
% dicom_sort(path1)
% Note: default path (UI_default_dir) can be edited for user's convenience.
%
% dicom_sort scrapes through dicom files in a single folder and create names:
% XXX(SCANNUMBER)_SCANNAME(WITHOUT UNDERSCORES OR GAPS!)
% Made for Siemens image conventions
% Raj Ramasawmy NHLBI, NIH, Dec 2017

UI_default_dir = '';

if nargin < 1
    path1 = uigetdir(UI_default_dir, 'Pick a folder containing dicoms');
end

%% Files need sorting
dir1 = dir(path1);

% Scrape through files & create names:
% XX(SCANNUMBER)_SCANNAME(WITHOUTUNDERSCORESORGAPS!)
disp('Step 1/3 - Scanning files');
for i = 3:length(dir1)
    
    % Check for Report Files
    temp = dir1(i).name;
    
    if regexp(temp(end-2:end), '.SR') == 1
        RepF(i) = 1;
    else
        RepF(i) = 0;
    end
    
    Name_ind{i} = 'blank'; % default
    
    % check if dicom file
    if ~isdir([path1 filesep dir1(i).name])
        if isdicom([path1 filesep dir1(i).name])
            % Extract Series Number and Name
            info1 = dicominfo([path1 filesep dir1(i).name]);
            Scan_ind(i) = info1.SeriesNumber;
            if(~isfield(info1, 'SeriesDescription'))
                info1.SeriesDescription = 'SCAN';
            end
            temp = regexprep(info1.SeriesDescription, '\W', ''); %Not Perfect!
            Name_ind{i} = [num2str(Scan_ind(i),'%02.f') '_' temp];
        end
    end
end

ScanVec = find(RepF == 0);
Name_ind{1} = 'blank';  Name_ind{2} = 'blank'; % ".", ".."
ScanNames = unique(Name_ind(ScanVec), 'stable');
ScanNames = ScanNames(2:end); % remove 'blank'
NumScans = length(ScanNames);

for i = 1:NumScans
    temp = ScanNames{i};
    over100(i) = ~isempty(str2num(temp(3)));
end

% Attempt to handle Phoenix report files <Siemens MRI specific?>
if(~isempty(find(RepF)))
    % If Report Files are included
    
    o100i = find(over100); % Make sure 99 (report files) is before 100
    if(~isempty(o100i))
        ScanNames2 = ScanNames;
        
        for i = 1:length(o100i)
            ScanNames{o100i(i)+1} = ScanNames2{o100i(i)};
        end
        
        ScanNames{o100i(1)} = '99_PhoenixReport';
        
        clear ScanNames2;
    else
        % if less than 100 scans
        ScanNames{NumScans+1} = '99_PhoenixReport';
    end
end
NumScans = length(ScanNames);

disp('Step 2/3 - Making directories');
% Make scan folders
for i = 1:length(ScanNames)
    mkdir(path1, ScanNames{i});
end

ScanNums = unique(Scan_ind(find(Scan_ind)) );

disp('Step 3/3 - Moving images');
% Sort Files
for i = 1:NumScans
    loc1 = find(Scan_ind == ScanNums(i));
    
    for j = 1:length(loc1)
        ScanFile = dir1(loc1(j)).name;
        movefile([path1 filesep ScanFile], [path1 filesep ScanNames{i}]);
        
    end
end

end