function [spikeanalysis] = InputOutput(ephysData, Cells, protocol_iteration)
%
% cells = number of cells you recorded in a single .dat bundle file. Usually it is one cell per .dat file
% protocol_iteration can be:
%   - omitted or []  → run ALL iterations that match protName
%   - a scalar (e.g. 2) → run just that iteration number among matches
%   - a vector (e.g. [1 3 5]) → run those specific iterations among matches
%
% For example:
%   [spikeanalysis] = InputOutput(ephysData, 1);          % all iterations
%   [spikeanalysis] = InputOutput(ephysData, 1, 1);       % only 1st iteration
%   [spikeanalysis] = InputOutput(ephysData, 1, [1 3 5]); % 1st, 3rd, 5th

%Sorting out traces and protocol info from .dat file for later analysis%

% ephysData is a structure array - a data type that allows you to group related data of different types 
% under one variable, with each field storing a separate piece of information
% A field in a structure is a named entity that stores data

allCells = fieldnames(ephysData); %Storing all field names in ephysData into a variable
disp('                       ');

% --- Preallocate output as empty and a running row index ---
spikeanalysis = [];
row = 0;

% --- store figures --- %
outDir = fullfile(pwd,'spiketrn_IFF_figs');
if ~exist(outDir,'dir'), mkdir(outDir); end

for iCell = Cells %Cells is our input argument, which is usually 1; However this number will increase if you record more than one cell in your .dat file

    cellName = allCells{iCell}; %Store name of the field that the cell is located in.

    %HARDCODED: Change to your protocol name
    % Look for pgf names similar to your protocol name and note their locations.
    protName = 'WRK_OSN'; 
    protLoc = find(strncmp(protName,ephysData.(cellName).protocols,length(protName))); % Store index of columns with name "WRK_OSN"

    % --- Choose which protocol indices to run ---
    if nargin < 3 || isempty(protocol_iteration)
        % Run all matching iterations (numbered relative to matches)
        iterNumbers = 1:numel(protLoc);                 % iteration numbers among matches
    else
        % Backward compatible: only the requested iteration(s)
        iterNumbers = protocol_iteration(:).';          % force row vector
    end

    % --- validate requested iteration numbers (relative to matches) ---
    if any(~isfinite(iterNumbers)) || any(iterNumbers~=floor(iterNumbers)) || any(iterNumbers<1)
        error('protocol_iteration must contain positive integer iteration numbers.');
    end
    if any(iterNumbers > numel(protLoc))
        error('Requested iteration(s) exceed available matches: requested max %d, available %d.', ...
              max(iterNumbers), numel(protLoc));
    end

    % Convert relative iteration numbers into absolute column indices
    protIdxList = protLoc(iterNumbers);                 % absolute indices into .data/.samplingFreq

    % --- loop over selected protocol iterations ---
    for k = 1:numel(protIdxList)
        selIdx  = protIdxList(k);       % absolute index into ephysData arrays
        iterTag = iterNumbers(k);       % relative iteration label (1..N among matches)

        selprotI = ephysData.(cellName).data{1, selIdx}; % all traces in this protocol iteration
        sampfreq = ephysData.(cellName).samplingFreq{1, selIdx}; % sampling frequency stored in .dat file

        total_datapoint = size(selprotI,1); % number of rows in ephys traces/number of datapoint = trace duration(s)*20,000
        total_sweep     = size(selprotI,2); % number of column in all ephys traces = number of traces

        close all

        pick = 1;     %%% set to 1 for manually setting (negative) threshold for spike detection
        scale = 10;   %%%% sets no of stds below mean for autothreshold

        prespikeno = 0;
        spikeno = 0;
        evoked_spike_ISIs = [];
        preISIs = [];

        % --- Per-iteration reset ---
        clear I t time

        %Looping through all the column (i.e. traces) - Analysis trace by trace%
        for sweep=1:total_sweep % nsweeps = number of trace/sweep in the current protocol
            bsl_duration = 6;
            I{sweep} = selprotI(:,sweep);   % Store all trace datapoint in the current trace into I{singlesweep}

            % Copying in Time data (normalising for sampling rate)
            for datapoint = 1:total_datapoint
                time(datapoint,sweep) = datapoint/sampfreq; %Create a 2D martix with row as time in 's' and column as current trace number
            end
            t{sweep} = time(:,sweep); % Select time for the current trace

            % === Top panel: raw trace with spikes ===
            fig = figure('Color','w');   % keep visible because you use ginput
            subplot(2,1,1);   % 2 rows, 1 column, 1st subplot (top)
            plot(t{sweep},I{sweep}); %plot time(s) in x-axis and current(A) in y-axis
            title(['WRK OSN Iter ' num2str(iterTag) ' (abs idx ' num2str(selIdx) '), Sweep ' num2str(sweep)]);
            minplot = min(I{sweep}); %#ok<NASGU>
            minaxis = min(I{sweep}) - 0.000000000005; %#ok<NASGU>
            axis([0 13.6 -0.0000000005 0.00000000005]); %HARDCODED:The limit of X-axis here is hardcoded.
            xlabel('Time (s)');
            ylabel('Current (A)');
            hold on

            % Selecting input threshold on the plot%
            if pick==1 
                disp('Click for spike threshold')
                [~,th] = ginput(1); % y-coordinate is the amplitude threshold
            else
                th = mean(alli)-scale.*std(alli); %#ok<NODEF>
                plot([allt(1) allt(end)],[mean(alli)-scale.*std(alli) mean(alli)-scale.*std(alli)],'m-') %#ok<NODEF>
            end

            % Ask for stimulation intensity
            stimulation_intensity = input('What is the stimulation intensity? ')*10; % convert to uA
            title(['WRK OSN Iter ' num2str(iterTag) ' (abs idx ' num2str(selIdx) '), stimulation ' num2str(stimulation_intensity) '𝛍A']);

            % pre-stim spikes: Background spikes%
            bsl_spike_timestamps = [];
            bsl_spike_indx = [];

            end_of_bsl_mask = t{sweep}<=bsl_duration;
            bsl_last_indx = find(end_of_bsl_mask, 1, 'last');

            for bsl_indx = 1:bsl_last_indx
                if I{sweep}(bsl_indx)>=th
                    if I{sweep}(bsl_indx+1)<th
                        bsl_spike_timestamps = [bsl_spike_timestamps; t{sweep}(bsl_indx)];
                        bsl_spike_indx = [bsl_spike_indx; bsl_indx];
                    end
                end
            end

            preTF = isempty(bsl_spike_timestamps); %#ok<NASGU>
            prespikeno = prespikeno + length(bsl_spike_timestamps);

            prespkt = [];
            prespka = [];
            prespki = [];

            %Looping through all the pre-stim spikes
            for m = 1:length(bsl_spike_indx)
                prespikew_i = bsl_spike_indx(m):bsl_spike_indx(m)+50; % ~2.5ms after crossing
                if prespikew_i(end)>(length(t{sweep}))
                    prespikew_i = bsl_spike_indx(m):(length(t{sweep}));  % for crosses right at the end of the sweep
                end

                prespikew_a = I{sweep}(prespikew_i);
                prespikew_t = t{sweep}(prespikew_i);
                [prespka(m), prespki(m)] = min(prespikew_a);
                prespkt(m) = prespikew_t(prespki(m));
            end

            % --- before filtering ---
            prespkt_filtered   = prespkt;      % default: nothing removed
            prespka_filtered   = prespka;
            prespkt_artefact   = [];
            prespka_artefact   = [];

            % --- filter spikes <1 ms apart ---
            if numel(prespkt) > 1
                q = find(diff(prespkt) < 0.001);
                if ~isempty(q)
                    arte_idx             = q + 1;
                    keep_idx             = setdiff(1:numel(prespkt), arte_idx);
                    prespkt_artefact     = prespkt(arte_idx);
                    prespka_artefact     = prespka(arte_idx);
                    prespkt_filtered     = prespkt(keep_idx);
                    prespka_filtered     = prespka(keep_idx);
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
                preISIs = diff(prespkt_filtered); %#ok<NASGU>
            end

            if ~isempty(prespkt_filtered)
                preTotal_SpikeNo = length(prespkt_filtered); %#ok<NASGU>

                if preTotal_SpikeNo > 1
                    prespk_trn_lth = prespkt_filtered(end) - prespkt_filtered(1); %#ok<NASGU>
                    if prespk_trn_lth > 0
                        preAve_Freq = preTotal_SpikeNo ./ prespk_trn_lth; % Hz
                    else
                        preAve_Freq = 0;
                    end

                elseif preTotal_SpikeNo == 1
                   prespk_trn_lth = 0; % Only 1 spike → cannot compute train length
                   preAve_Freq = preTotal_SpikeNo ./ bsl_duration;
                end
            else
                preTotal_SpikeNo = 0; %#ok<NASGU>
                prespk_trn_lth = 0; %#ok<NASGU>
                preAve_Freq    = 0;
            end

            %post-stim spikes: response spikes%
            evoked_spike_timestamps = [];
            evoked_spike_indx = [];

            beginning_of_recording_mask = t{sweep}>(bsl_duration+0.01); % Start of response window; +10 ms to avoid detecting stimulation artefact
            recording_start_indx = find(beginning_of_recording_mask, 1, 'first');

            for recording_indx = recording_start_indx:length(I{sweep})-1
                if I{sweep}(recording_indx)>=th
                    if I{sweep}(recording_indx+1)<th % -1 to prevent indexing outside of the vector
                        evoked_spike_timestamps = [evoked_spike_timestamps; t{sweep}(recording_indx)];
                        evoked_spike_indx = [evoked_spike_indx; recording_indx];
                    end
                end
            end

            TF = isempty(evoked_spike_timestamps); %#ok<NASGU>
            spikeno = spikeno + length(evoked_spike_timestamps); %#ok<NASGU>

            spkt = [];
            spka = [];
            spki = [];

            %Looping through all the evoked spikes
            for m = 1:length(evoked_spike_indx)
                spikew_i = evoked_spike_indx(m):evoked_spike_indx(m)+50; % ~2.5ms after crossing
                if spikew_i(end)>(length(t{sweep}))
                    spikew_i = evoked_spike_indx(m):(length(t{sweep}));
                end
                spikew_a = I{sweep}(spikew_i);
                spikew_t = t{sweep}(spikew_i);
                [spka(m), spki(m)] = min(spikew_a);
                spkt(m) = spikew_t(spki(m));
            end

            % Update graph axis base on pre-spikes or evoked-spike amplitude
            if ~isempty(spkt)
                ymin = min(spka) - 20e-12; % 20 pA
                ymax = min(spka) + 100e-12; % 100 pA
                axis([0 13.6 ymin ymax]);
            elseif isempty(spkt) && ~isempty(prespka)
                ymin = min(prespka) - 20e-12; % 20 pA
                ymax = min(prespka) + 100e-12; % 100 pA
                axis([0 13.6 ymin ymax]);
            else
                axis([0 13.6 -0.0000000005 0.00000000005]);
            end

            % --- before filtering ---
            spkt_filtered   = spkt;      % default: nothing removed
            spka_filtered   = spka;
            spkt_artefact   = [];
            spka_artefact   = [];

            % --- filter spikes <1 ms apart ---
            if numel(spkt) > 1
                artefact_contaminated_spike_indx = find(diff(spkt)<.001);
                if ~isempty(artefact_contaminated_spike_indx)
                    evoked_spike_artefact_indx = artefact_contaminated_spike_indx + 1; % spike_artefact_indx = spike_indx + 1
                    evoked_spike_keep_indx = setdiff(1:numel(spkt),evoked_spike_artefact_indx); % Remove artefact indexes 
                    spkt_filtered   = spkt(evoked_spike_keep_indx);
                    spka_filtered   = spka(evoked_spike_keep_indx);
                    spkt_artefact   = spkt(evoked_spike_artefact_indx);
                    spka_artefact   = spka(evoked_spike_artefact_indx);
                end
            end

            % Calculate IFF and determine spike train based on IFF
            evoked_spike_ISIs_vector = diff(spkt_filtered);
            evoked_IFF_vector = 1./evoked_spike_ISIs_vector;

            % check if all IFF is =< bsl average firing frequency
            if ~isempty(evoked_IFF_vector) && all(preAve_Freq >= evoked_IFF_vector)
                % --- If true, no spikes are included
                spkt_trn = [];
                spka_trn = [];
            else
                % --- Find first indx when evoked_IFF_vector <= preAve_Freq
                idxDrop = find(evoked_IFF_vector <= preAve_Freq, 1, 'first');
                % --- If none, include all spikes
                if isempty(idxDrop)
                    spkt_trn = spkt_filtered;
                    spka_trn = spka_filtered;
                else
                    % include spikes before IFF FIRST dropped/equal to preAve_Freq
                    spkt_trn = spkt_filtered(1:idxDrop);
                    spka_trn = spka_filtered(1:idxDrop);
                end
            end

            TF2 = isempty(spkt_trn); %#ok<NASGU>

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
                evoked_spike_ISIs = diff(spkt_trn); %#ok<NASGU>
            end
            hold off

            % === Bottom panel: IFF ===
            subplot(2,1,2);   % 2 rows, 1 column, 2nd subplot (bottom)
            if length(spkt_trn) > 1
                plot(spkt_filtered(2:end), evoked_IFF_vector, 'o-');
                hold on
                if exist('idxDrop','var') && ~isempty(idxDrop)
                    xline(spkt_filtered(idxDrop+1),'r--'); % Beginning (inclusive) of spike exclusion
                end
                yline(preAve_Freq,'r--');
                hold off
                xlabel('Time (s)');
                ylabel('IFF (Hz)');
                title('Instantaneous Firing Frequency');
            else
                text(0.5,0.5,'Not enough spikes for calculating IFF','Units','normalized',...
                    'HorizontalAlignment','center');
            end

            if ~isempty(evoked_spike_timestamps) && ~isempty(spkt_trn)
                Total_SpikeNo = length(spkt_trn);
                Total_Time = t{sweep}(end) - (bsl_duration+0.01);
                spkt_trn_lth = spkt_trn(end) - spkt_trn(1);
                if spkt_trn_lth>0
                    Ave_Freq = length(spkt_trn) ./ spkt_trn_lth;
                else
                    Ave_Freq = 0;
                end
                if numel(spkt_trn)>1
                    tmpISI = diff(spkt_trn);
                    CV = std(tmpISI) ./ mean(tmpISI);
                    Max_Freq = 1 / min(tmpISI);
                else
                    CV = NaN;
                    Max_Freq = 0;
                end
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

            % Subtracting evoked spikes by preAve_Freq*spkt_trn_lth
            bsl_spike_in_spkt_trn  = round(preAve_Freq*spkt_trn_lth);
            bsl_subtracted_SpikeNo = Total_SpikeNo - bsl_spike_in_spkt_trn;
            if bsl_subtracted_SpikeNo < 0
                bsl_subtracted_SpikeNo = 0;
            end

            % Subtracted spike section
            row = row + 1;
            spikeanalysis(row, 1)  = iterTag;                  % relative iteration number among matches
            spikeanalysis(row, 2)  = stimulation_intensity;    % uA
            spikeanalysis(row, 3)  = sweep;                    % sweep index
            spikeanalysis(row, 4)  = bsl_subtracted_SpikeNo;   % baseline-subtracted spikes in train
            spikeanalysis(row, 5)  = Total_SpikeNo;            % raw spikes in train
            spikeanalysis(row, 6)  = preAve_Freq;              % baseline avg freq (Hz)
            spikeanalysis(row, 7)  = spkt_trn_lth;             % train length (s)
            spikeanalysis(row, 8)  = bsl_spike_in_spkt_trn;    % expected baseline spikes in train

            % Post-stimulus spike section
            spikeanalysis(row, 9)  = sweep;                    % (duplicate of col 3)
            spikeanalysis(row,10)  = Total_SpikeNo;            % total spikes
            spikeanalysis(row,11)  = Total_Time;               % analyzed duration (s)
            spikeanalysis(row,12)  = spkt_trn_lth;             % train length (s)
            spikeanalysis(row,13)  = Ave_Freq;                 % avg frequency (Hz)
            spikeanalysis(row,14)  = Max_Freq;                 % max frequency (Hz)
            spikeanalysis(row,15)  = CV;                       % coefficient of variation of ISIs
            spikeanalysis(row,16)  = Latency;                  % latency (s)

            % Pre-stimuli spike section
            if exist('preTotal_SpikeNo','var') && exist('prespk_trn_lth','var') && exist('preAve_Freq','var')
                spikeanalysis(row,17)  = preTotal_SpikeNo;        
                spikeanalysis(row,18)  = prespk_trn_lth;           
                spikeanalysis(row,19)  = preAve_Freq;              
            else
                spikeanalysis(row,17:19) = [0 0 0];
            end

            disp("press any key to continue")

            % === SAVE FIGURE / SUBPLOTS ===
            drawnow;  % make sure graphics are up-to-date
            fnameBase = sprintf('WRKOSN %02d_Sweep%02d_%g_uA', iterTag, sweep, stimulation_intensity);
            pngPath = fullfile(outDir, [fnameBase '.png']);
            exportgraphics(fig, pngPath, 'Resolution', 300);
            pause;
            close(fig);
        end % sweep
    end % protocol iteration
end % cell

% Group by iteration (col 1), then sort by stimulation intensity (col 2)
if ~isempty(spikeanalysis)
    spikeanalysis = sortrows(spikeanalysis, [1 2]);
end
