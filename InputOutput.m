function [spikeanalysis] = InputOutput(ephysData, Cells, protocol_iteration)

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

    cellName = allCells{iCell}; %Store name of the field that the cell is located in.

    %HARDCODED: Change to your protocol name
    % Look for pgf names similar to your protocol name and note their locations.
    protName = 'WRK_OSN'; 
    protLoc = find(strncmp(protName,ephysData.(cellName).protocols,length(protName))); % Store index of columns with name "WRK_OSN"
                                                                                       
    % Look for pgf names similar to your protocol name and note their locations.
    
    %protocol_iteration = protocol; %protocol is our input arguement, which is to select the number of protocol iteration you want to analyse
    
    selprotI = ephysData.(cellName).data{1,protLoc(protocol_iteration)}; % assessing row 1 and column protLoc(j)
                                                       % i.e. all ephys traces in that protocol iteration 

    total_datapoint = size(selprotI,1); % number of rows in ephys traces/number of datapoint = trace duration(s)*20,000
    total_sweep = size(selprotI,2); % number of column in all ephys traces = number of traces
    sampfreq = ephysData.(cellName).samplingFreq{1,protLoc(protocol_iteration)}; % sampling frequency stored in .dat file
    
    close all
   
    pick = 1;  %%% set to 1 for manually setting (negative) threshold for spike detection
    scale = 10;  %%%% sets no of stds below mean for autothreshold
    
    prespikeno = 0;
    spikeno = 0;
    evoked_spike_ISIs = [];
    preISIs = [];
    
end

%Looping through all the column (i.e. traces) - Analysis trace by trace%

for sweep=1:total_sweep % nsweeps = number of trace/sweep in the current protocol
    
    bsl_duration = 6;
    I{sweep} = selprotI(:,sweep);   % Store all trace datapoint in the current trace into I{singlesweep}
    
    % Copying in Time data (normalising for sampling rate)
    for datapoint = 1:total_datapoint
        time(datapoint,sweep) = datapoint/sampfreq; %Create a 2D martix with row as time in 's' and column as current trace number
                                                    %E.g. 260020th datapoint/20K Hz = 13.001 sec
    end
    
    t{sweep} = time(:,sweep); % Select time for the current trace
    
    plot(t{sweep},I{sweep}); %plot time(s) in x-axis and current(A) in y-axis
    title(['Sweep ' num2str(sweep)]);
    minplot = min(I{sweep});
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

    bsl_spike_timestamps = []; % Create an empty list for storing timepoint (s) of spikes
    bsl_spike_indx = [];  % Create an empty list for index of spikes

    end_of_bsl_mask = t{sweep}<=bsl_duration; % Logical output: "1" for all time in t{sweep} if it is before 6s.
                                  
    bsl_last_indx = find(end_of_bsl_mask, 1, 'last'); % Find the last index of the baseline window that is non-zero

    for bsl_indx = 1:bsl_last_indx
        if I{sweep}(bsl_indx)>=th %%% so if current trace datapoint in baseline window is above thresh...
            if I{sweep}(bsl_indx+1)<th %%% ...but next point is below thresh... This is because spike has a downward deflection
                
                bsl_spike_timestamps = [bsl_spike_timestamps; t{sweep}(bsl_indx)];% Timestamping spike and append in pre_spikecrosses
                                                               % vertical concatenation. This means the code 
                                                               % takes the existing array pre_spikecrosses 
                                                               % (which initially is empty) and appends the 
                                                               % new time value as a new row at the end

                bsl_spike_indx = [bsl_spike_indx; bsl_indx];        % Same as above but for index
            
            end
        end
    end
    
    preTF = isempty(bsl_spike_timestamps); %Logical output: whether pre_spikecrosses is empty or not
    
    prespikeno = prespikeno + length(bsl_spike_timestamps); % Store the number of spikes in baseline window
    
    prespkt = []; % Empty list for storing all pre-spikes time(s)
    prespka = []; % for storing lowest peak amplitude of all pre-spikes
    prespki = []; % for storing the index of all pre-spikes
    
    %Looping through all the pre-stim spikes
    for m = 1:length(bsl_spike_indx)
        prespikew_i = bsl_spike_indx(m):bsl_spike_indx(m)+50; % HARDCODED: spike window ~2.5ms after crossing
                                                                % Store index for the entire spike width
                                                                % Qs: Why 2.5 ms
        if prespikew_i(end)>(length(t{sweep}))
            prespikew_i = bsl_spike_indx(m):(length(t{sweep}));  % for crosses right at the end of the sweep
        end
        
        prespikew_a = I{sweep}(prespikew_i);           % Store current (A) of the entire spike trough - from the point of spike threshold to 2.5 ms later
        prespikew_t = t{sweep}(prespikew_i);           % ... for time of the entire spike trough
        [prespka(m), prespki(m)] = min(prespikew_a);   % prespka = storing the lowest current (A) value of the spike trough
                                                       % prespki = storing the index of lowest current (A) value of the spike trough
        prespkt(m) = prespikew_t(prespki(m));          % prespkt = storing the time(s) of the lowest current (A) value of the spike trough
        
    end

    % --- before filtering ---
    prespkt_filtered   = prespkt;      % default: nothing removed
    prespka_filtered   = prespka;
    prespkt_artefact   = [];           % default: no artefacts; reset every sweeps
    prespka_artefact   = [];
    
    % --- filter spikes <1 ms apart ---
    if numel(prespkt) > 1

        % Find spike indx where they are less that 1ms apart
        q = find(diff(prespkt) < 0.001);
        if ~isempty(q)
            arte_idx             = q + 1;                 % the "artefact" spikes
                                                          % q+1: List of index of pre-stim spike that's <1ms away from their previous spike
            keep_idx             = setdiff(1:numel(prespkt), arte_idx);   % Remove indexes in q+1 from 1:length(prespkt)
                                                                          % Remove index of artefact spikes from pre-stim spikes
    
            prespkt_artefact     = prespkt(arte_idx);     % Time of spike artefacts
            prespka_artefact     = prespka(arte_idx);     % Amplitude of spike artefacts
            prespkt_filtered     = prespkt(keep_idx);     % Time of filtered spikes
            prespka_filtered     = prespka(keep_idx);     % Amplitude of filtered spikes
        end
    end
    
    % --- plot baseline spike ---
    if ~isempty(prespkt_filtered)
        plot(prespkt_filtered, prespka_filtered, 'ro');   % pre-stim spikes
    end
    
    % --- plot baseline spike artefact ---
    if ~isempty(prespkt_artefact)
        plot(prespkt_artefact, prespka_artefact, 'bo');   % artefacts
    end
    
    % --- storing baseline spike ISI ---
    if length(prespkt_filtered)>1
        preISIs = diff(prespkt_filtered);
    end

    if ~isempty(prespkt_filtered)
    preTotal_SpikeNo = length(prespkt_filtered)

        if preTotal_SpikeNo > 1
            prespk_trn_lth = prespkt_filtered(end) - prespkt_filtered(1)
    
            if prespk_trn_lth > 0
                preAve_Freq = preTotal_SpikeNo ./ prespk_trn_lth % Hz
            else
                preAve_Freq = 0
            end
        
    
        elseif preTotal_SpikeNo == 1
           prespk_trn_lth = 0 % Only 1 spike → cannot compute train length
           preAve_Freq = preTotal_SpikeNo ./ bsl_duration
        
        end
        

    else
            preTotal_SpikeNo = 0
            prespk_trn_lth = 0
            preAve_Freq    = 0
    end

    % %Filter spikes artefact that's 1ms apart from each other
    % if length(prespkt)>1
    %     q = find(diff(prespkt)<.001);  % q+1: List of index of pre-stim spike that's <1ms away from their previous spike
    %     if ~isempty(q)
    %         cdi = setdiff([1:length(prespkt)],q+1);  % Remove indexes in q+1 from [1:length(prespkt)]
    %                                                  %%% so removing any 'spikes' within 1ms of each other from the list of prespkt
    % 
    %         prespkt_filtered = prespkt(cdi); % 
    %         prespka_filtered = prespka(cdi);
    %         prespkt_artefact = prespkt(q+1);
    %         prespka_artefact = prespka(q+1);
    %     end
    % end
    % 
    % % Circling out pre-stim spikes (filtered)
    % if ~isempty(prespkt_filtered)
    %     plot(prespkt_filtered,prespka_filtered,'ro')
    %     plot(prespkt_artefact,prespka_artefact,'bo','MarkerFaceColor','b')
    % 
    %     if length(prespkt_filtered)>1
    %         preISIs = diff(prespkt_filtered);
    %     end
    % end

    
%post-stim spikes: response spikes%

    evoked_spike_timestamps = []; % Create an empty list for storing timepoint (s) of spike
    evoked_spike_indx = []; % ... for storing spike index
    

    beginning_of_recording_mask = t{sweep}>(bsl_duration+0.01); % Start of response window; +10 ms to avoid detecting stimulation artefact

    recording_start_indx = find(beginning_of_recording_mask, 1, 'first'); % find the index of start time point of recording window

    for recording_indx = recording_start_indx:length(I{sweep})-1
        if I{sweep}(recording_indx)>=th % so if point is above thresh...
            if I{sweep}(recording_indx+1)<th % ...but next point is below thresh...
                                             % -1 to prevent indexing outside of the vector
                evoked_spike_timestamps = [evoked_spike_timestamps; t{sweep}(recording_indx)];        % Timestamping spike and append in evoked_spike_timestamps
                                                               % vertical concatenation. This means the code 
                                                               % takes the existing array spikecrosses 
                                                               % (which initially is empty) and appends the 
                                                               % new time value as a new row at the end

                evoked_spike_indx = [evoked_spike_indx; recording_indx];                % Same as above but for index
            end
        end
    end
    
    
   TF = isempty(evoked_spike_timestamps); % Logical output "0": evoked_spike_timestamps is not empty
    
   spikeno = spikeno + length(evoked_spike_timestamps); % Count the number of spikes
    
   spkt = [];  % Empty list for storing spike time(s)
   spka = [];  % for storing lowest peak amplitude of spike
   spki = [];  % for storing the index of spike
    
    %Looping through all the evoked spikes
    for m = 1:length(evoked_spike_indx)
        spikew_i = evoked_spike_indx(m):evoked_spike_indx(m)+50; %%% spike window ~2.5ms after crossing
        if spikew_i(end)>(length(t{sweep}))
            spikew_i = evoked_spike_indx(m):(length(t{sweep}));  % for crosses right at the end of the sweep
        end
        
        spikew_a = I{sweep}(spikew_i);            % Store current (A) of the entire spike trough - from the point of spike threshold to 2.5 ms later
        spikew_t = t{sweep}(spikew_i);            % ... for time of the entire spike trough
        [spka(m), spki(m)] = min(spikew_a);       % spka = storing the lowest current (A) value of the spike trough
        spkt(m) = spikew_t(spki(m));              % spki = ...index of...
                                                  % spkt = ...time(s) of...

        
    end

    % Update graph axis base on pre-spikes or evoked-spike amplitude
    if ~isempty(spkt)
        ymin = min(spka) - 20e-12; % 20 pA
        ymax = min(spka) + 100e-12; % 40 pA
        axis([0 13.6 ymin ymax]);

    elseif isempty(spkt) & ~isempty(prespka)
        ymin = min(prespka) - 20e-12; % 20 pA
        ymax = min(prespka) + 100e-12; % 40 pA
        axis([0 13.6 ymin ymax]);

    else
        axis([0 13.6 -0.0000000005 0.00000000005]); %HARDCODED:The limit of X-axis here is hardcoded.
    end
    

     % --- before filtering ---
    spkt_filtered   = spkt;      % default: nothing removed
    spka_filtered   = spka;
    spkt_artefact   = [];        % default: no artefacts; reset every sweeps
    spka_artefact   = [];
    
    
     % --- filter spikes <1 ms apart --- There will be spike artefact, if the spike threshold are set too close resting current.
    if numel(spkt) > 1
        artefact_contaminated_spike_indx = find(diff(spkt)<.001);  %%% so looking for spike crosses less than 1ms apart...
        if ~isempty(artefact_contaminated_spike_indx)
            evoked_spike_artefact_indx = artefact_contaminated_spike_indx + 1; % spike_artefact_indx = spike_indx + 1

            evoked_spike_keep_indx = setdiff(1:numel(spkt),evoked_spike_artefact_indx); % Remove artefact indexes 

            spkt_filtered   = spkt_filtered(evoked_spike_keep_indx); % Time of filtered spike
            spka_filtered   = spka(evoked_spike_keep_indx);          % Amplitude(A) of filtered spike
            spkt_artefact   = spkt(evoked_spike_artefact_indx);      % Time of spike artefact
            spka_artefact   = spka(evoked_spike_artefact_indx);      % Amplitude of spike artefact

        end
    end

    % Calculate IFF and determine spike train based on IFF

    evoked_spike_ISIs_vector = diff(spkt_filtered);
    evoked_IFF_vector = 1./evoked_spike_ISIs_vector;

     % check if all IFF is =< bsl average firing frequency
    if all(preAve_Freq >= evoked_IFF_vector)
    
        % --- If true, no spikes are included
        spkt_trn = [];
        spka_trn = [];

    else
        % --- Find first indx when evoked_IFF_vector <= preAve_Freq
        idxDrop = find(evoked_IFF_vector <= preAve_Freq, 1, 'first');

        % --- If none, include all spikes
        if isempty(idxDrop)
            spkt_trn = spkt_filtered;
            spka_trn = spka_filtered
    
        % --- else, include spikes before IFF FIRST dropped/equal to preAve_Freq
        % --- spike with IFF == preAve_Freq is excluded
        else
            spkt_trn = spkt_filtered(1:idxDrop);
            spka_trn = spka_filtered(1:idxDrop)
        end
    
        %numel(spkt_beforeDrop)

    end
    
    % disp('Click on end of spike train')
    % 
    % [x1,~] = ginput(1);  % Store the x-coordinate(i.e. time(s)) of your input
    % intrain_mask = spkt_filtered<x1; % Logical output of whether spike are within your selected time(s)
    % 
    % spkt_trn = spkt_filtered(intrain_mask); % Time(s) of spike within your selected time
    % spka_trn = spka_filtered(intrain_mask)           % Amplitude(A)....
    % 
    
    TF2 = isempty(spkt_trn);

     % --- green dots for evoked spike ---
     if ~isempty(spkt_trn)
         plot(spkt_trn,spka_trn,'go')
     end
    % --- purple dots for evoked spike artefact ---
     if ~isempty(spkt_artefact)
         plot(spkt_artefact,spka_artefact,'o','Color', [0.5 0 0.5])
     end
    % --- store evoked spike ISI ---
     if length(spkt_trn)>1
         evoked_spike_ISIs = diff(spkt_trn);
     end

     hold off
    
    %figure(2)
    %hist(ISIs,50)
    %disp(mean(ISIs))
    %disp(std(ISIs))
    
     
    if TF == 0 && TF2 == 0
        Total_SpikeNo = length(spkt_trn);
        Total_Time = t{sweep}(end) - (bsl_duration+0.01);
        spkt_trn_lth = spkt_trn(end) - spkt_trn(1);
        Ave_Freq = length(spkt_trn) ./ spkt_trn_lth;
        CV = std(evoked_spike_ISIs) ./ mean(evoked_spike_ISIs);
        Max_Freq = 1 / min(evoked_spike_ISIs);
        Latency = evoked_spike_timestamps(1) - (bsl_duration+0.01);
        
    else
        Total_SpikeNo = 0;
        Total_Time = t{sweep}(end) - (bsl_duration+0.01);
        spkt_trn_lth = 0;
        Ave_Freq = 0;
        CV = NaN;
        Max_Freq = 0;
        Latency = NaN;
    end
    
    spikeanalysis(sweep, 1) = sweep;
    spikeanalysis(sweep, 2) = Total_SpikeNo;
    spikeanalysis(sweep, 3) = Total_Time;
    spikeanalysis(sweep, 4) = spkt_trn_lth;
    spikeanalysis(sweep, 5) = Ave_Freq;
    spikeanalysis(sweep, 6) = Max_Freq;
    spikeanalysis(sweep, 7) = CV;
    spikeanalysis(sweep, 8) = Latency;
    
    spikeanalysis(sweep, 9) = preTotal_SpikeNo;
    spikeanalysis(sweep, 10) = prespk_trn_lth;
    spikeanalysis(sweep, 11) = preAve_Freq;
       
    
    disp("press any key to continue")
    pause;

end

