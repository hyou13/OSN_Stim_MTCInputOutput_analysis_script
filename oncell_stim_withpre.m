function [spikeanalysis] = oncell_stim_withpre(ephysData, Cells, protocol)

% cells = number of cells you recorded in a single .dat bundle file. Usually it is one cell per .dat file
% protocol = number of repetition.
% For example, if you repeated the same protocol for three times, then you need to run the function thrice
% E.g. First repeat of cell 1 - (ephysData, 1,1);2nd repeat - (ephysData, 1,2);3rd repeat - (ephysData, 1,3)

%Sorting out traces and protocol info from .dat file for later analysis%

% ephysData is a structure array - a data type that allows you to group related data of different types 
% under one variable, with each field storing a separate piece of information
% A field in a structure is a named entity that stores data

allCells = fieldnames(ephysData); %Storing all field names in ephysData into a variable
disp('                       ');


for iCell = Cells %Cells is our input argument, which is usually 1; However this number will increase if you record more than one cell in your .dat file

    cellName = allCells{iCell}; %storing cell into a variable

    %HARDCODED: Change to your protocol name
    % Look for pgf names similar to your protocol name and note their locations.
    protName = 'WRK_OSN stim'; 
    protLoc = find(strncmp(protName,ephysData.(cellName).protocols,length(protName))); % Find index of columns with name "WRK_OSN"
                                                                                       %Qs: Why do we need bracket here?
    % Look for pgf names similar to your protocol name and note their locations.
    
    j = protocol; %protocol is our input arguement, which is to select the protocol you want to analyse
    
    selprotI = ephysData.(cellName).data{1,protLoc(j)}; % assessing row 1 and column protLoc(j)
                                                       % i.e. all ephys traces in the protcol

    s = size(selprotI,1); % number of rows in ephys traces/number of datapoint = trace duration(s)*20,000
    nsweeps = size(selprotI,2); % number of column in all ephys traces = number of traces
    sampfreq = ephysData.(cellName).samplingFreq{1,protLoc(j)}; % sampling frequency stored in .dat file
    
    close all
   
    pick = 1;  %%% set to 1 for manually setting (negative) threshold for spike detection
    scale = 10;  %%%% sets no of stds below mean for autothreshold
    
    prespikeno = 0;
    spikeno = 0;
    ISIs = [];
    preISIs = [];
    
end

%Looping through all the column (i.e. traces) - Analysis trace by trace%

for k=1:nsweeps % nsweeps = number of trace/sweep in the current protocol
    
    I{k} = selprotI(:,k);   % Store all trace datapoint in the current trace into I{k}
    
    % Copying in Time data (normalising for sampling rate)
    for n = 1:s
        time(n,k) = n/sampfreq; %Create a 2D martix with row as time in 's' and column as current trace number
    end
    
    t{k} = time(:,k); % Select time for the current trace
    
    plot(t{k},I{k}); %plot time(s) in x-axis and current(A) in y-axis
    title(['Sweep ' num2str(k)]);
    minplot = min(I{k});
    minaxis = minplot - 0.000000000005;
    axis([0 13.6 -0.0000000005 0.00000000005]); %HARDCODED:The limit of X-axis here is hardcoded.
    
    hold on
    
%Selecting input threshold on the plot%

    if pick==1 
        disp('Click for spike threshold')
        [~,th] = ginput(1); % ginput(n) - identify coordinates of 'n' point and store them into variable [x,y]
                            % y-coordinate is store as 'thr' which is the amplitude threshold for spike ; x-coordinate is neglected
    else
        th = mean(alli)-scale.*std(alli);
        %plot([allt(1) allt(end)],[mean(alli) mean(alli)],'m-')
        plot([allt(1) allt(end)],[mean(alli)-scale.*std(alli) mean(alli)-scale.*std(alli)],'m-')
    end
        
%pre-stim spikes: Background spikes%

    pre_spikecrosses = []; % Create an empty list for storing timepoint (s) of spikes
    pre_spikecrossi = [];  % Create an empty list for index of spikes

    b = t{k}<6; % HARDCODED: At 6s, we gave our electrical stimulation.
                   % Essentially this line stores all the time(s) of our baseline window into a variable.

    index_endpt = find(b, 1, 'last'); % Find the last index of the baseline window

    for j = 1:index_endpt
        if I{k}(j)>=th %%% so if trace datapoint in baseline window is above thresh...
            if I{k}(j+1)<th %%% ...but next point is below thresh... This is because spike has a downward deflection
                
                pre_spikecrosses = [pre_spikecrosses; t{k}(j)];% Timestamping spike and append in pre_spikecrosses
                                                               % vertical concatenation. This means the code 
                                                               % takes the existing array pre_spikecrosses 
                                                               % (which initially is empty) and appends the 
                                                               % new time value as a new row at the end

                pre_spikecrossi = [pre_spikecrossi; j];        % Same as above but for index
            
            end
        end
    end
    
    preTF = isempty(pre_spikecrosses); %Logical output: whether pre_spikecrosses is empty or not
    
    prespikeno = prespikeno + length(pre_spikecrosses); % Store the number of spikes in baseline window
    
    prespkt = []; % Empty list for storing pre-spike time(s)
    prespka = []; % for storing lowest peak amplitude of pre-spike
    prespki = []; % for storing the index of pre-spike
    
    %Looping through all the pre-stim spikes
    for m = 1:length(pre_spikecrossi)
        prespikew_i = pre_spikecrossi(m):pre_spikecrossi(m)+50; % HARDCODED: spike window ~2.5ms after crossing
                                                                % Store index for the entire spike width
                                                                % Qs: Why 2.5 ms
        if prespikew_i(end)>(length(t{k}))
            prespikew_i = pre_spikecrossi(m):(length(t{k}));  % for crosses right at the end of the sweep
        end
        
        prespikew_a = I{k}(prespikew_i); % Store current (A) of the entire width of spike
        prespikew_t = t{k}(prespikew_i); % for time
        [prespka(m), prespki(m)] = min(prespikew_a);   %storing the current (A) value of the peak of the spike
        prespkt(m) = prespikew_t(prespki(m)); %storing the time(s) of the peak of the spike
        
    end
    
    if length(prespkt)>1
        q = find(diff(prespkt)<.001);  %%% so looking for spike crosses less than 1ms apart...
        if ~isempty(q)
            cdi = setdiff([1:length(prespkt)],q+1); %%% so removing any 'spikes' within 1ms of each other. Qs: Why do we need to do this?
            prespkt = prespkt(cdi); 
            prespka = prespka(cdi);
        end
    end
    
    % Qs: Never saw this hapenning
    if ~isempty(prespkt)
        plot(prespkt,prespka,'ro')
        if length(prespkt)>1
            preISIs = diff(prespkt);
        end
    end
    
    
%post-stim spikes: response spikes%

    spikecrosses = []; % Create an empty list for storing timepoint (s) of spike
    spikecrossi = []; % ... for storing spike index

    a = t{k}>6.01; % HARDCODED: Start of response window

    index_startpt = find(a, 1, 'first'); % find the index of start time point of recording window

    for j = index_startpt:length(I{k})-1 %Qs: why -1?
        if I{k}(j)>=th %%% so if point is above thresh...
            if I{k}(j+1)<th %%% ...but next point is below thresh...
                spikecrosses = [spikecrosses; t{k}(j)];        % Timestamping spike and append in spikecrosses
                                                               % vertical concatenation. This means the code 
                                                               % takes the existing array spikecrosses 
                                                               % (which initially is empty) and appends the 
                                                               % new time value as a new row at the end

                spikecrossi = [spikecrossi; j];                % Same as above but for index
            end
        end
    end
    
    
   TF = isempty(spikecrosses); % Logical output: whether pre_spikecrosses is empty or not
    
    spikeno = spikeno + length(spikecrosses); % Count the number of spikes
    
    spkt = [];  % Empty list for storing spike time(s)
    spka = [];  % for storing lowest peak amplitude of spike
    spki = [];  % for storing the index of spike
    
    for m = 1:length(spikecrossi)
        spikew_i = spikecrossi(m):spikecrossi(m)+50; %%% spike window ~3ms after crossing
        if spikew_i(end)>(length(t{k}))
            spikew_i = spikecrossi(m):(length(t{k}));  % for crosses right at the end of the sweep
        end
        
        spikew_a = I{k}(spikew_i);            % Store current (A) of the entire width of spike
        spikew_t = t{k}(spikew_i);            % ... for time ...
        [spka(m), spki(m)] = min(spikew_a);   % storing the current (A) value of the peak of the spike
        spkt(m) = spikew_t(spki(m));          % storing the time(s) of the peak of the spike. Qs: If spki(m) is a current(a) value, how can it be used to index out the time?

        
    end
    
    
    if length(spkt)>1
        q = find(diff(spkt)<.001);  %%% so looking for spike crosses less than 2ms apart...
        if ~isempty(q)
            cdi = setdiff([1:length(spkt)],q+1);
            spkt = spkt(cdi); %%% so removing any 'spikes' within 2ms of each other
            spka = spka(cdi);
        end
    end
    
%stopped here%
    
disp('Click on end of spike train')

[x1,~] = ginput(1);  %%% plot rough location of last peak
intrain = spkt<x1;

spkt_trn = spkt(intrain);
spka_trn = spka(intrain)

TF2 = isempty(spkt_trn);

 if ~isempty(spkt_trn)
        plot(spkt_trn,spka_trn,'go')
        if length(spkt_trn)>1
            ISIs = diff(spkt_trn);
        end
 end
    
    hold off
    
    %figure(2)
    %hist(ISIs,50)
    %disp(mean(ISIs))
    %disp(std(ISIs))
    
    
    if preTF == 0
        preTotal_SpikeNo = length(prespkt);
        prespk_trn_lth = prespkt(end) - prespkt(1);
        preAve_Freq = preTotal_SpikeNo ./ prespk_trn_lth;

    else
        preTotal_SpikeNo = 0;
        prespk_trn_lth = 0;
        preAve_Freq = 0;
        
    end
    
     
    if TF == 0 && TF2 == 0
        Total_SpikeNo = length(spkt_trn);
        Total_Time = t{k}(end) - 6.01; %HARDCODED
        spkt_trn_lth = spkt_trn(end) - spkt_trn(1);
        Ave_Freq = length(spkt_trn) ./ spkt_trn_lth;
        CV = std(ISIs) ./ mean(ISIs);
        Max_Freq = 1 / min(ISIs);
        Latency = spikecrosses(1) - 6.01; %HARDCODED
        
    else
        Total_SpikeNo = 0;
        Total_Time = t{k}(end) - 6.01; %HARDCODED
        spkt_trn_lth = 0;
        Ave_Freq = 0;
        CV = NaN;
        Max_Freq = 0;
        Latency = NaN;
    end
    
    spikeanalysis(k, 1) = k;
    spikeanalysis(k, 2) = Total_SpikeNo;
    spikeanalysis(k, 3) = Total_Time;
    spikeanalysis(k, 4) = spkt_trn_lth;
    spikeanalysis(k, 5) = Ave_Freq;
    spikeanalysis(k, 6) = Max_Freq;
    spikeanalysis(k, 7) = CV;
    spikeanalysis(k, 8) = Latency;
    
    spikeanalysis(k, 9) = preTotal_SpikeNo;
    spikeanalysis(k, 10) = prespk_trn_lth;
    spikeanalysis(k, 11) = preAve_Freq;
       
    
    disp("press any key to continue")
    pause;
end


