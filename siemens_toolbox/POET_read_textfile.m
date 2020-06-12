function [Data] =POET_read_textfile(xdata)
% [Data] =POET_read_textfile(xdata <optional>)
% 
% Can print out any combination of RF, ADC, Gradients from POET
% Data ordering will have to be sorted by user (header is printed out)

if nargin < 1
   [a,b] = uigetfile('*.*');
   xdata = [b filesep a]; clear a b;
end

%% Initialise Data

f = fopen(xdata);
stop = 0;
while stop==0;
    s = fgetl(f);
    if strfind(s,'Start Time')>0
        i = strfind(s,'=');
        i2 = strfind(s,'SRT');
        starttime = str2num(s(i+1:i2-1));
    end
    
    if strfind(s,'End Time')>0
        i = strfind(s,'=');
        i2 = strfind(s,'SRT');
        endtime = str2num(s(i+1:i2-1));
    end
    
    if strfind(s,'Time Step')>0
        i = strfind(s,'=');
        i2 = strfind(s,'SRT');
        timestep = str2num(s(i+1:i2-1));
    end
    
    if( ~isempty(strfind(s,'Signal')) || ~isempty(strfind(s,'Gradient')))
        a = textscan(s, '%s', 'Delimiter', '\t');
        numChannels = length(a{1});
        a{1} % print to screen the channel header and Data order
        stop =1;
    end
    
end

fclose(f);
numPoints = ((endtime - starttime)/timestep);
Data = zeros(numPoints, numChannels);

%% Recognize labels
POET_data_names = {'ADC','RF', 'X ', 'Y ', 'Z '}; 
% Header names (typical)
%     'ADC Signal Data'
%     'RF-Signal  (ch. 0, 1H, ##.# MHz)'
%     'X Gradient (GPA 0)'
%     'Y Gradient (GPA 0)'
%     'Z Gradient (GPA 0)'
%     'Numeric Crystal Oscillator'
channel_order = zeros(1,numChannels);
for i = 1:numChannels
    channel_name = a{1}{i};
    for j = 1:length(POET_data_names)
        temp = regexp(channel_name, POET_data_names{j});
        if ~isempty(temp)
            channel_order(j) = i;
        end
    end
end
channel_order(find(channel_order==0))= [];

%% Grab Data

f = fopen(xdata);
stop = 0; n =1;
while stop==0;
    s = fgetl(f);
    if ( ~isempty(strfind(s,'Signal')) || ~isempty(strfind(s,'Gradient')))
        s = fgetl(f); % Empty Row after Header
        stop2 =0;
        % Begin grabbing
        while stop2 ==0;
            s = fgetl(f);
            if s==-1;
                stop2 =1; % end of file
            else
                b = textscan(s, '%s', 'Delimiter', '\t');
                for j = 1:length(b{1})
                    temp = cell2mat(b{1}(j));
                    if(~isempty(temp)) % Textscan factors empty columns
                        Data(n,j) = str2num(temp);
                    end
                end
            end
            n = n+1;
        end
        stop =1;
    end
end

fclose(f);

%% Chop off excess Points (for debugging edited txt files)
% if n < numPoints
%    Data = Data(1:n,:);
% end

%% Plot Channels (optional)

% === plot ADC,RF,X,Y,Z ===
% make_dev;
% dev.plot_poet_data(Data);

% === plot all ===
for i = 1:numChannels
    channelNames{i} = a{1}{numChannels + 1- i}(isstrprop(a{1}{numChannels + 1- i}, 'alpha')); % reverse order to match plot order in POET.
end
scale_pd = max(abs(Data),[],1);

yax_des = 1:2.1:(1 + 2.1*numChannels-1);
figure, hold on;
for i = 1:numChannels
    plot(Data(:,numChannels + 1- i)/scale_pd(numChannels + 1- i) + yax_des(i), 'k-', 'LineWidth', 1.5)
end
ytick('off'); ylim([-0.1 max(yax_des)+1.1]);
set(gca,'ytick',yax_des,'yticklabel',channelNames)
xlabel('Time (us)');

plot_invert = 0;
if plot_invert
    % plot for presentation with black background:
    %  set(gcf, 'OuterPosition', [150 500 1500 800])
    
    lines = findobj(gcf,'Type','Line');
    for i = 1:numel(lines)
        lines(i).LineWidth = 1.5; % may want to thicken more
        lines(i).Color = 'w';
    end
    set(gca, 'Color', 'k')
    a = gcf;
    a.InvertHardcopy = 'off';
    
end

end





