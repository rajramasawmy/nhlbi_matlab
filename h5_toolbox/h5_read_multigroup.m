function [outputStruct] = h5_read_multigroup(filename, dataSelect)
% Load all data
% [a] = h5_read_multigroup('csm_transfer\out.h5');
% Select data to load
% [b] = h5_read_multigroup('csm_transfer\out.h5', 0);
% Load this data (i.e. > 0)
% [c] = h5_read_multigroup('csm_transfer\out.h5', 9);

%%
if nargin < 1
    [a,b] = uigetfile('*.*');
    filename = [b filesep a]; clear a b;
end
a = h5info(filename);

numRecons = length(a.Groups);
if nargin > 1
    if dataSelect == 0
        disp(sprintf('************\nChoose Data \n************'));
        for nRc = 1:numRecons
            b = a.Groups(nRc).Groups;
            disp(['Image set: ' num2str(nRc) '   Recon ID: ' a.Groups(nRc).Name])
        end
        disp('************');
        dataSelect = input('Choose Image Set: ');
        %     dataSelect = str2double(dataSelect);
    end
    loopVector = dataSelect;
else
    loopVector = 1:numRecons;
end

%%
data_list = {};
for nRc = loopVector
    b = a.Groups(nRc).Groups;
    disp(['Loading Image set: ' num2str(nRc) '   Recon ID: ' a.Groups(nRc).Name])
    data_list{nRc} = a.Groups(nRc).Name;
    
    for i = 1 : length(b)
        temp = b(i).Name;
        image_name = regexp(temp, 'image_');
        image_name = temp(image_name:end);
        
        data = squeeze(double(h5read(filename, [temp '/data'])));
        
        if i == 1 % export header once
            info = h5read(filename, [temp '/header']);
        end
        
        outputStruct.(matlab.lang.makeValidName(['data' num2str(nRc)])).(matlab.lang.makeValidName(image_name)) = data;
        
    end
    outputStruct.(matlab.lang.makeValidName(['data' num2str(nRc)])).header = info;
end

outputStruct.data_list = data_list;

end