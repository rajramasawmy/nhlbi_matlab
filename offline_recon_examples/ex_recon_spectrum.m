%%-------------------------------%%
%%----MRI reconstruction code----%%
%%-------------------------------%%
% function [freq, stats] = recon_spectrum(file)
% Siemens service RF noise 
% input:    file = ISMRMRD .h5 data file 
%
% output:   freq = all spectra
    %       stats = per-channel summary

% R Ramasawmy May 2019 NHLBI 


function [freq, stats, header] = recon_spectrum(file)
%% Read data file
make_nhlbi_toolbox;

file = nhlbi_toolbox.run_path_on_sys(file); % incorporate with NHLBI toolbox

raw_data= h5read(file, '/dataset/data');
ismrmrd_s = read_h5_header(file); disp(' ');disp('### Protocol Name ###');disp(ismrmrd_s.measurementInformation.protocolName);disp(' ');
header = ismrmrd_s;

samples = double(raw_data.head.number_of_samples(1));
channels = double(raw_data.head.active_channels(1));

pe1 = double(max(raw_data.head.idx.kspace_encode_step_1))+1;
pe2 = double(max(raw_data.head.idx.kspace_encode_step_2))+1;
averages = double(max(raw_data.head.idx.average))+1;
slices = double(max(raw_data.head.idx.slice))+1;
contrasts = double(max(raw_data.head.idx.contrast))+1;
phases = double(max(raw_data.head.idx.phase))+1;
reps = double(max(raw_data.head.idx.repetition)+1);
sets = double(max(raw_data.head.idx.set))+1;
dt = raw_data.head.sample_time_us(1)*1e-6;

disp(' ');disp('### Experiment Dimensions ###');disp(' ');
Experiment_parameters = {'Samples', 'PE1', 'PE2', 'Averages', 'Slices', 'Contrasts', 'Phases', 'Repetitions', 'Sets', 'Channels'}';
Value = [samples pe1 pe2 averages slices contrasts phases reps sets channels]';
disp(table( Experiment_parameters,Value )); clear Experiment_parameters Value; disp(' ');

%% List coil names

asi_names = fieldnames(ismrmrd_s.acquisitionSystemInformation);

solid_int = 0;
coil_label = cell(channels,2);

for i = 1:length(asi_names)
    if regexp(asi_names{i}, 'coilLabel')
        solid_int = solid_int + 1;
        coil_label{solid_int,1} = ismrmrd_s.acquisitionSystemInformation.(matlab.lang.makeValidName(asi_names{i})).coilNumber;
        coil_label{solid_int,2} = ismrmrd_s.acquisitionSystemInformation.(matlab.lang.makeValidName(asi_names{i})).coilName;
    end
end

Data_CoilNum = coil_label(:,1);
Data_CoilName = coil_label(:,2);
disp(table(Data_CoilNum, Data_CoilName)); disp(' ');

clear Data_CoilName Data_CoilNum 

%% Grab data

kspace = complex(zeros([samples pe1 pe2 averages slices contrasts phases reps sets channels],'single'));
% disp(['Kspace dims: ' num2str(size(kspace))])

for ii = 1:length(raw_data.data)
    
    d1 = raw_data.data{ii};
    d2 = complex(d1(1:2:end), d1(2:2:end));
    d3 = reshape(d2, samples, channels); %  RE & IM (2)
    
    kspace(:,...
        raw_data.head.idx.kspace_encode_step_1(ii)+1, ...
        raw_data.head.idx.kspace_encode_step_2(ii)+1, ...
        raw_data.head.idx.average(ii)+1, ....
        raw_data.head.idx.slice(ii)+1 , ...
        raw_data.head.idx.contrast(ii)+1, ...
        raw_data.head.idx.phase(ii)+1, ...
        raw_data.head.idx.repetition(ii)+1, ...
        raw_data.head.idx.set(ii)+1, ...
        :) = d3;
    
end

%% FFT 
% Currently spits out all coil data (typically ran with body coil)
% PCA and choose a single, primary coil?

kspace = squeeze(kspace);
kspace = reshape(kspace, [samples reps*averages channels]);

freq = abs(ismrm_transform_kspace_to_image(kspace,1));

Hz_limit = 0.5/dt;
kHz_limit = Hz_limit/1e3;
kHz_range = -kHz_limit:0.5*kHz_limit:kHz_limit;

%% Plot

sp1 = floor(sqrt(channels));
sp2 = ceil(channels/sp1);

figure,
for i = 1:channels
    subplot(sp1,sp2, i); 
    shadedErrorBar(1:samples, mean(freq(:,:,i),2)', std(freq(:,:,i),[],2)'); title(coil_label(i,2), 'Interpreter', 'none'); 
    
    % 2016b onwards
    %     xticks(1+round(linspace(0,samples-1,5)))
    %     xticklabels(kHz_range); xlabel('kHz'); title(['Channel ' num2str(i)]);
    
    stats(i,:) = [mean(mean(freq(:,:,i),2)) mean( std(freq(:,:,i),[],2) )]; % dont think this is informative.. 
end

%%

header.coil_label = coil_label(:,2);

end

