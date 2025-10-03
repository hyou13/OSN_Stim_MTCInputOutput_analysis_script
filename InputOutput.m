function [spikeanalysis] = InputOutput(ephysData, Cells, protocol_iteration)

% === Output ===
% 1. summaryTbl of the cell - Treatment_wpi_ExptCond_Date_CellNumber_Expter (e.g. MMZ_2w_aCSF-GBZ-NBQX_APV_250409_cell2_HM) 
% 2. spikeanalysis - raw data of the cell
% 3. LinearRegression_slope_plots
% 4. spiketrn_IFF_figs


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

% === Define output diretory ===
AllOutput_Dir = '/Users/hyou/Documents/OSN_Stim_MTCInputOutput_analysis_script/Output_240424to0525';

allCells = fieldnames(ephysData); %Storing all field names in ephysData into a variable
disp('                       ');

% --- Preallocate output as empty and a running row index ---
spikeanalysis = [];
row = 0;

% --- store figures --- %
IFF_outDir = fullfile(AllOutput_Dir,'spiketrn_IFF_figs_240424to0525');
if ~exist(IFF_outDir,'dir'), mkdir(IFF_outDir); end

for iCell = Cells %Cells is our input argument, which is usually 1; However this number will increase if you record more than one cell in your .dat file

    % === ask for metadata and build the output filename ===
    Treatment    = strtrim(input('Enter Treatment (e.g., PBS or MMZ): ', 's'));
    wpi          = strtrim(input('Enter week post injection (e.g., 2w): ', 's'));
    date         = strtrim(input('Enter experiment date (e.g., yy/mm/dd 250423): ', 's'));  
    cell_number  = strtrim(input('Enter cell number (e.g., cell3): ', 's'));
    experimenter = strtrim(input('Enter experimenter name (e.g. HM, Lorcan): ', 's')); 

    % create file for individual cell IFF_plot
    IFF_filename = sprintf('%s_%s_%s_%s_%s', ...
                     Treatment, wpi, date, cell_number, experimenter);
    IFF_OutDir = fullfile(IFF_outDir, IFF_filename);
    if ~exist(IFF_OutDir,'dir')
        mkdir(IFF_OutDir);
    end

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

            %disp("press any key to continue")

            % === SAVE FIGURE / SUBPLOTS ===
            drawnow;  % make sure graphics are up-to-date
            fnameBase = sprintf('WRKOSN %02d_Sweep%02d_%g_uA', iterTag, sweep, stimulation_intensity);
            pngPath = fullfile(IFF_OutDir, [fnameBase '.png']);
            exportgraphics(fig, pngPath, 'Resolution', 300);
            %pause;
            close(fig);
        end % sweep
    end % protocol iteration
end % cell

% Group by iteration (col 1), then sort by stimulation intensity (col 2)
if ~isempty(spikeanalysis)
    spikeanalysis = sortrows(spikeanalysis, [1 2]);
end

% ===== RUNDOWN ANALYSIS =====

% all distinct intensities, except 3000 & 10000
exclude = [3000 10000];
intensity_exclusion_mask = ~ismember(spikeanalysis(:,2), exclude);
intensities = unique(spikeanalysis(intensity_exclusion_mask, 2));                  
nI          = numel(intensities);

% --- For intensity with repeated sweep in each iterations, calculate sweep-averaged spike no ---
perI_SpikeVals   = cell(nI,1);
iter_counts      = zeros(nI,1);

% Loop through all stimulation intensities
for k = 1:nI
    % extract all spikeanalysis rows for this intensity
    thisI = intensities(k);
    rowsK = spikeanalysis(spikeanalysis(:,2) == thisI, :);

    % sort by iteration then sweep
    rowsK = sortrows(rowsK, [1 3]);
    
    iterIdx   = rowsK(:,1);     % list of iteration (contains repeat)
    spikeVal  = rowsK(:,4);     % bsl_subtracted_spikeno

    % list of unique iteration (without repeats)
    iterList = unique(iterIdx, 'stable');

    % create table for storing sweep-averaged bsl_subtracted_spikeno of each iteration
    avg_spikeVals = NaN(numel(iterList),1);

    % --- average spike across sweeps within the same iteration --- 
    for ii = 1:numel(iterList)
        % create a mask for each iteration number
        repeated_iteration_mask       = (iterIdx == iterList(ii));
        % averaging across sweeps
        avg_spikeVals(ii)= round(mean(spikeVal(repeated_iteration_mask), 'omitnan'));
    end

    % store sweep-averaged spike for each stimulation intensity
    perI_SpikeVals{k}  = avg_spikeVals;               % vector used to build SpikesTbl
    
    % vector of iteration list for each stimulation intensity
    perI_iters{k} = iterList;              % keep the iteration order for mapping rows

    %Number of iteration for each stimulation intensity
    iter_counts(k)    = size(perI_iters{k},1);
end

% === Table A: sweep-average bsl_subtracted_SpikeNo (of each ordered iteration) x intensities ===
% missing values are filled with NaN
allIters = unique(spikeanalysis(:,1), 'sorted');   % all iteration IDs present (e.g., [1 2 3 4 5])
nRows    = numel(allIters);
SpikesTbl = NaN(nRows, nI);                        % rows = iterations (without repeats), cols = intensities

% --- map sweep-average bsl_subtracted_SpikeNo to each iteration in SpikesTbl ---
for k = 1:nI
    if ~isempty(perI_iters{k})
        % tf(i)  = true if perI_iters{k} exists in allIters, else false
        % ps(i)  = row index within allIters where iters_k(i) appears; 0 if not found
        [tf, pos] = ismember(perI_iters{k}, allIters);

        % place averaged values in the correct iteration rows
        SpikesTbl(pos(tf), k) = perI_SpikeVals{k}(tf);
        % pos(tf)           = matching iteration between perI_iters and allIters, rows in SpikesTbl
        % k                 = stimulation intensity, columns in SpikesTbl
        % perI_vals{k}(tf)  = iteration-macthing sweep-average bsl_subtracted_SpikeNo
    
    end

end

% === Table B: Percentage of rundown ===
RunDownTbl = NaN(size(SpikesTbl));

for k = 1:nI
    % for the current stimulation intensity...
    spike_per_iteration_vector = SpikesTbl(:, k);
    % identify first row that is non-NaN
    validrow_idx = find(~isnan(spike_per_iteration_vector));           

    % pick the first NON-ZERO among the valid rows 
    first_non0idx = find(spike_per_iteration_vector(validrow_idx) ~= 0, 1, 'first');   % position within validrow_idx
        
        if isempty(first_non0idx)
            % all valid values are zero → no usable baseline; skip this column
            continue
        end
    
    % baseline as first non-NaN and non-0 row
    baseline = spike_per_iteration_vector(validrow_idx(first_non0idx));

    % divide all valid rows by the baseline (exclude earlier NaN and 0 rows)
    later = validrow_idx(first_non0idx:end);   % keep all non-NaN rows
    RunDownTbl(later, k) = spike_per_iteration_vector(later) ./ baseline;

end


% --- caculate iteration-averaged rundown to both table ---
avgRundown_per_iteration = mean(RunDownTbl(:,1:nI), 2, 'omitnan');
RunDownTbl(:, nI+1) = avgRundown_per_iteration;
SpikesTbl(:, nI+1) = avgRundown_per_iteration;


% === Sectioning SpikesTbl into different experimental condition === 

% number of rows in SpikesTbl/PctTbl
nRows_SpikesTbl = size(SpikesTbl,1);

% % === ask for metadata and build the output filename ===
% Treatment    = strtrim(input('Enter Treatment (e.g., PBS or MMZ): ', 's'));
% wpi          = strtrim(input('Enter week post injection (e.g., 2w): ', 's'));
% date         = strtrim(input('Enter experiment date (e.g., yy/mm/dd 250423): ', 's'));  
% cell_number  = strtrim(input('Enter cell number (e.g., cell3): ', 's'));
% experimenter = strtrim(input('Enter experimenter name (e.g. HM, Lorcan): ', 's')); 

% --- Ask user for experimental condition names (in order) ---
conditions = strtrim(input(['Enter experimental conditions in order, ' ...
                     'separated by commas or spaces (e.g., aCSF, GBZ, NBQX_APV): '], 's'));

if isempty(conditions)
    error('No condition names entered.');
end

% Split on commas/semicolons or whitespace, then trim
if contains(conditions, ',') || contains(conditions, ';')
    parts = regexp(conditions, '\s*[;,]\s*', 'split');
else
    parts = strsplit(conditions);   % split on whitespace
end

% drop empty cells
parts = strtrim(parts);
parts = parts(~cellfun('isempty', parts));

% (Optional) de-duplicate but keep first occurrence
[~, ia] = unique(lower(parts), 'stable');
parts   = parts(ia);

% --- Build the cell array of condition names ---
condNames_array = parts;   % e.g., {'aCSF','GBZ','NBQX_APV'}
cond_cutoffs    = zeros(1,numel(condNames_array));   % inclusive end-row per condition (0 = skip)
prev_condEnd    = 0;

% Join the selected condition names with hyphens (keep underscores within names)
if iscell(condNames_array)
    cond_string = strjoin(condNames_array, '-');
else
    cond_string = string(condNames_array);
end
cond_string = char(cond_string);


for i = 1:numel(condNames_array)
    condEnd = NaN;
    % ask until valid
    % DO NOT accept: NaN,smaller than preEnd, larger than nRows & integer
    while isnan(condEnd) || condEnd < prev_condEnd || condEnd > nRows_SpikesTbl || condEnd ~= floor(condEnd)
        % formatted string
        condEnd = input(sprintf('Up until which iteration is %s? (>= %d & <= %d; enter %d to skip): ', ...
                             condNames_array{i}, prev_condEnd+1, nRows_SpikesTbl, prev_condEnd));
    end
    % 0 or >= prevEnd
    cond_cutoffs(i) = condEnd;     
    % Update prev_condEnd for next condition
    prev_condEnd    = condEnd;
end

% Convert cond_cutoffs → row-index ranges (cells)
condRows = cell(1,numel(condNames_array));
prev_condEnd  = 0;

for i = 1:numel(condNames_array)
    if cond_cutoffs(i) > prev_condEnd
        % Take all the number between prev_condEnd+1 to cuttoff number
        % aCSF cutoff = 2; take 0+1 to 2 → [1 2]
        condRows{i} = (prev_condEnd+1):cond_cutoffs(i);
    else
        condRows{i} = [];  % skipped
    end
    prev_condEnd = cond_cutoffs(i);
end

% Pack into a struct for easy access (e.g. condIdx.aCSF)
condIdx = cell2struct(condRows, condNames_array, 2);

% === column average spikes within conditions if rundown is >0.75 ====

% --- Ask user for experimental condition names (in order) ---
colAvg_conditions = strtrim(input( ...
    ['(Rundown dependent) Enter experimental conditions to be '...
     'calculated for column average ' ...
     'in order, separated by commas or spaces ' ...
     '(e.g., aCSF, GBZ, NBQX_APV): '], 's'));

if isempty(colAvg_conditions)
    error('No condition names entered.');
end

% Split on commas/semicolons or whitespace, then trim
if contains(colAvg_conditions, ',') || contains(colAvg_conditions, ';')
    colAvg_parts = regexp(colAvg_conditions, '\s*[;,]\s*', 'split');
else
    colAvg_parts = strsplit(colAvg_conditions);   % split on whitespace
end

% drop empty cells
colAvg_parts = strtrim(colAvg_parts);
colAvg_parts = colAvg_parts(~cellfun('isempty', colAvg_parts));

% (Optional) de-duplicate but keep first occurrence
[~, ia] = unique(lower(colAvg_parts), 'stable');
colAvg_parts   = colAvg_parts(ia);

% --- Build the cell array of condition names ---
colAvg_condNames_array = colAvg_parts;   % e.g., {'aCSF','GBZ','NBQX_APV'}

% specify column index for spikes and rundown
SpikeColsindx        =   nI;
AvgRundownColsindx   =   nI + 1;

IO_averaged = struct();              % will store per-condition 1×nI row vectors
avgRundown_by_condition = struct();  % (optional) keep the gating metric too
IO_drop = struct();

% create figure for linear regression of IO slope
fig = figure;
hold on;
legendHandles = [];
legendEntries = {};

% Create output folder for storing plots (consistent name)
LinearRegression_slope_plots_outDir = fullfile(AllOutput_Dir, 'LinearRegression_slope_plots_240424to0525');
if ~exist(LinearRegression_slope_plots_outDir,'dir'), mkdir(LinearRegression_slope_plots_outDir); end

for i = 1:numel(colAvg_condNames_array)
    
    % storing name & rows for each condition
    name = colAvg_condNames_array{i};
    rows = condIdx.(name);
    
    % if rows is empty skip to next condition
    if isempty(rows), continue; end %% Need to add an empty array

    % --- Get per-iteration avgRundown for this condition and keep those >0.75 ---
    % column of %rundown for these rows
    condRundown = SpikesTbl(rows, AvgRundownColsindx);

    % row passes if rundown > 0.75
    keepMask    = ~isnan(condRundown) & condRundown > 0.75;
    dropMask    = isnan(condRundown) | condRundown < 0.75;

    rows_kept = rows(keepMask);
    rows_drop = rows(dropMask);
    
    % if no rows other than baseline row passes the criterion for this condition...
    if isempty(rows_kept) || numel(rows_kept)<=1
        % put column average is empty
        avgRundown_by_condition.(name) = NaN;
        IO_averaged.(name)             = [];

        % --- store drop rows ---
        IO_drop.(name) = SpikesTbl(rows_drop,:);

        % skips everything below and continue to the next condition
        fprintf('%s: Not enough row (0 or <=1) passed the > 0.75 criterion.\n', name);
        continue;
    end
    
    % --- average rundown for condition (column-wise) ---
    cond_AvgRundown = mean(SpikesTbl(rows_kept, AvgRundownColsindx), 'omitnan');
    
    % --- average spikes across kept rows only (column-wise) ---
    cond_colAvg = round(mean(SpikesTbl(rows_kept, 1:SpikeColsindx), 1, 'omitnan'));

    % --- Stimulation threshold: find first non-NaN/0 index --- % 

    % mask of entries that are not NaN and not zero
    mask = ~isnan(cond_colAvg) & (cond_colAvg ~= 0);
    % first index that passes the mask
    first_indx = find(mask, 1, 'first');
    input_threshold = intensities(first_indx);

    % --- Omax: as the mean of 0.8-1.0 of maximum output ----
    Omax_mask = (cond_colAvg./max(cond_colAvg))>=0.8;
    Omax = round(mean(cond_colAvg(Omax_mask)));
    
    % --- slope value & plotting (overlaid for all conditions) ---
    % Inputs (vectors)
    X    = log10(intensities(:));   % log10 of input (Nx1)
    Ymat = cond_colAvg(:);          % outputs (N x 1)
    slope = NaN;                    % default in case we can't fit
    
    % set up a color per condition so data+fit share the same color
    ax = gca; 
    colorOrder = get(ax, 'ColorOrder');
    nColors = size(colorOrder,1);
    
    % keep paired non-NaN samples
    keep = ~isnan(Ymat) & ~isnan(X);
    Xf = X(keep);
    Yf = Ymat(keep);

    % Sort by x so lines connect left-to-right
    [x_sorted, sortIdx] = sort(Xf);
    y_sorted = Yf(sortIdx);

    % ----- Find the last zero (or closest-to-zero) output on the sorted data -----
    zeroIdx = find(y_sorted == 0, 1, 'last');
    if isempty(zeroIdx)
        [~, zeroIdx] = min(abs(y_sorted));
    end

    fitRange = zeroIdx:numel(y_sorted);
    x_fit = x_sorted(fitRange);
    y_fit = y_sorted(fitRange);

    % color for this condition
    thisColor = colorOrder(mod(i-1, nColors) + 1, :);

    % Plot raw data for this condition
    hData = plot(x_sorted, y_sorted, 'o-', 'Color', thisColor, 'MarkerFaceColor', thisColor);
    legendHandles(end+1) = hData; %#ok<AGROW>
    legendEntries{end+1} = ['Data-' name]; %#ok<AGROW>
    
    % Fit & plot linear regression (if enough points)
    if numel(x_fit) >= 2
        p = polyfit(x_fit, y_fit, 1);   % y = p(1)*x + p(2)
        slope = p(1);

        x_line = linspace(min(x_fit), max(x_fit), 100);
        y_line = polyval(p, x_line);
        hFit = plot(x_line, y_line, '--', 'LineWidth', 1.5, 'Color', thisColor);

        legendHandles(end+1) = hFit; %#ok<AGROW>
        legendEntries{end+1} = ['linear fit-' name]; %#ok<AGROW>
        legend(legendHandles, legendEntries, 'Location', 'best');
    else
        fprintf('%s: Not enough points in fit window to compute linear fit.\n', name);
    end

    IO_averaged.(name) = [cond_colAvg, cond_AvgRundown, input_threshold, Omax, slope];
    IO_drop.(name)     = SpikesTbl(rows_drop,:);

end

% Finalize figure
title('Linear regression of I-O plot (log10 transformed): All conditions');
xlabel('log10(Intensity)'); ylabel('spikes');
hold off;

% Save a single overlaid plot
filename = sprintf('%s_%s_%s_%s_%s_LinearRegression_overlaid.png', ...
                 Treatment, wpi, date, cell_number, experimenter);

fpath = fullfile(LinearRegression_slope_plots_outDir, filename);
exportgraphics(fig, fpath, 'Resolution', 300);
close(fig);

% === column average spikes within conditions without considering rundown (e.g. NBQX_APV) ====

% --- Ask user for experimental condition names (in order) ---
Rdwn_indep_colAvg_conditions = strtrim(input( ...
    ['(Rundown independent) Enter experimental conditions to be calculated ' ...
     'for column average in order, separated by ' ...
     'commas or spaces (e.g., NBQX_APV). ' ...
     'To skip, press enter: '], 's'));

Rdwn_indep_IO_averaged = struct();   % will store per-condition 1×nI row vectors

if isempty(Rdwn_indep_colAvg_conditions)
    Rdwn_indep_cond_AvgRundown = NaN;
    Rdwn_indep_cond_colAvg = NaN(1, SpikeColsindx);
    Rdwn_indep_IO_averaged.NoSelection = [Rdwn_indep_cond_colAvg,Rdwn_indep_cond_AvgRundown];

else

    % Split on commas/semicolons or whitespace, then trim
    if contains(Rdwn_indep_colAvg_conditions, ',') || contains(Rdwn_indep_colAvg_conditions, ';')
        Rdwn_indep_colAvg_parts = regexp(Rdwn_indep_colAvg_conditions, '\s*[;,]\s*', 'split');
    else
        Rdwn_indep_colAvg_parts = strsplit(Rdwn_indep_colAvg_conditions);   % split on whitespace
    end
    
    % drop empty cells
    Rdwn_indep_colAvg_parts = strtrim(Rdwn_indep_colAvg_parts);
    Rdwn_indep_colAvg_parts = Rdwn_indep_colAvg_parts(~cellfun('isempty', Rdwn_indep_colAvg_parts));
    
    % --- Build the cell array of condition names ---
    Rdwn_indep_colAvg_condNames_array = Rdwn_indep_colAvg_parts;   % e.g., {'aCSF','GBZ','NBQX_APV'}
    
    
    for i = 1:numel(Rdwn_indep_colAvg_condNames_array)
        
        % storing name & rows for each condition
        Rdwn_indep_name = Rdwn_indep_colAvg_condNames_array{i};
        Rdwn_indep_rows = condIdx.(Rdwn_indep_name);
        
        % if rows is empty skip to next condition
        if isempty(Rdwn_indep_rows), continue; end
        
        % --- average spikes across kept rows only (column-wise) ---
        Rdwn_indep_cond_AvgRundown = mean(SpikesTbl(Rdwn_indep_rows, AvgRundownColsindx), 'omitnan');
        Rdwn_indep_cond_colAvg = round(mean(SpikesTbl(Rdwn_indep_rows, 1:SpikeColsindx), 1, 'omitnan'));
        Rdwn_indep_IO_averaged.(Rdwn_indep_name) = [Rdwn_indep_cond_colAvg,Rdwn_indep_cond_AvgRundown];
        
    end
end

% ===== Build the summary table and export to Excel =====

% 0) (Optional) show zeros as blanks in Excel
Sp = SpikesTbl;      Sp(Sp==0) = NaN;
Rd = RunDownTbl;     Rd(Rd==0) = NaN;

% 1) Column names = intensities + avgRundown_per_iteration (+ optional extras)
nI        = size(SpikesTbl,2) - 1;                   % last col is avgRundown
baseNames = [compose('I_%g', intensities(:).'), "avgRundown_per_iteration"];
baseNames = cellstr(matlab.lang.makeValidName(baseNames));

% If IO_averaged vectors are wider than SpikesTbl (e.g., Threshold/Rmax/pchipSlope),
% add extra column names automatically (or rename below as you like).
ioNames     = fieldnames(IO_averaged);
extraWidth  = 0;
for i = 1:numel(ioNames)
    extraWidth = max(extraWidth, numel(IO_averaged.(ioNames{i})) - numel(baseNames));
end
extraWidth  = max(extraWidth, 0);
if extraWidth > 0
    defaultExtra = ["Threshold","Rmax","linear_interpolation_Slope"];                   % customize if needed
    if numel(defaultExtra) < extraWidth
        defaultExtra = [defaultExtra, compose("Extra_%d", 1:(extraWidth-numel(defaultExtra)))];
    end
    extraNames = cellstr(defaultExtra(1:extraWidth));
else
    extraNames = {};
end

varNames  = [baseNames, extraNames];
nColsOut  = numel(varNames);

% helpers
padrow = @(r) [r(:).', NaN(1, max(0, nColsOut-numel(r)))];              % pad/truncate to width
labels = {};                                                            % first column (text)
nums   = NaN(0, nColsOut);                                              % numeric block we append to
nanRow = NaN(1, nColsOut);

% 2) IO averaged section
labels{end+1,1} = 'IO averaged'; nums(end+1,:) = nanRow;
for i = 1:numel(ioNames)
    name = ioNames{i};
    row  = IO_averaged.(name);
    labels{end+1,1} = name;      nums(end+1,:) = padrow(row);
end
labels{end+1,1} = '';            nums(end+1,:) = nanRow;

% 3) Rundown-independent IO averaged section
io2Names = fieldnames(Rdwn_indep_IO_averaged);
labels{end+1,1} = 'Rdwn_indep_IO_averaged'; nums(end+1,:) = nanRow;
for i = 1:numel(io2Names)
    name = io2Names{i};
    row  = Rdwn_indep_IO_averaged.(name);
    labels{end+1,1} = name;      nums(end+1,:) = padrow(row);
end
labels{end+1,1} = '';            nums(end+1,:) = nanRow;

% 4) IO dropped section (just list sizes)
dropNames = fieldnames(IO_drop);
labels{end+1,1} = 'IO dropped';  nums(end+1,:) = nanRow;
for i = 1:numel(dropNames)
    name   = dropNames{i};
    [r,c]  = size(IO_drop.(name));
    labels{end+1,1} = sprintf('%s: [%dx%d double]', name, r, c);
    nums(end+1,:)   = nanRow;                                        % keep numeric area blank
end

% 5) SpikesTbl section
labels{end+1,1} = 'SpikesTbl';  nums(end+1,:) = nanRow;
for r = 1:size(Sp,1)
    labels{end+1,1} = '';        nums(end+1,:) = padrow(Sp(r,:));
end
labels{end+1,1} = '';            nums(end+1,:) = nanRow;                 % blank separator

% 6) RunDownTbl section
labels{end+1,1} = 'RunDownTbl'; nums(end+1,:) = nanRow;
for r = 1:size(Rd,1)
    labels{end+1,1} = '';        nums(end+1,:) = padrow(Rd(r,:));
end
labels{end+1,1} = '';            nums(end+1,:) = nanRow;

% 7) Build table and write to Excel
summaryTbl = array2table(nums, 'VariableNames', varNames);
summaryTbl = addvars(summaryTbl, labels, 'Before', 1, 'NewVariableNames', 'Section');

% Build filename: Treatment_wpi_cond1-cond2-..._date_cellX_(Experimenter).xlsx
% Create output folder for storing summaryTbl (consistent name)
summaryTbl_outDir = fullfile(AllOutput_Dir, 'summaryTbl_240424to0525');
if ~exist(summaryTbl_outDir,'dir'), mkdir(summaryTbl_outDir); end

outFile = sprintf('%s_%s_%s_%s_%s_%s.xlsx', ...
                  Treatment, wpi, cond_string, date, cell_number, experimenter);
outXlsx = fullfile(summaryTbl_outDir, outFile);

writetable(summaryTbl, outXlsx, 'FileType','spreadsheet');
fprintf('Wrote %s\n', outXlsx);
% Output spikesanalysis file

% Create output folder for storing summaryTbl (consistent name)
spikeanalysis_outDir = fullfile(AllOutput_Dir, 'spikeanalysis_240424to0525');
if ~exist(spikeanalysis_outDir,'dir'), mkdir(spikeanalysis_outDir); end

spikeanalysis_outFile = sprintf('%s_%s_%s_%s_%s_%s_spikeanalysis.xlsx', ...
                  Treatment, wpi, cond_string, date, cell_number, experimenter);

spikeanalysis_outXlsx = fullfile(spikeanalysis_outDir,spikeanalysis_outFile);

% Column headers for the 19 columns in spikeanalysis (must match order above)
colHeaders = { ...
    'iterTag', 'stimulation_intensity', 'sweep', 'bsl_subtracted_SpikeNo', ...
    'Total_SpikeNo', 'preAve_Freq', 'spkt_trn_lth', 'bsl_spike_in_spkt_trn', ...
    'post_sweep', 'post_Total_SpikeNo', 'post_Total_Time', 'post_spkt_trn_lth', ...
    'post_Ave_Freq', 'post_Max_Freq', 'post_CV', 'Latency', ...
    'preTotal_SpikeNo', 'prespk_trn_lth', 'preAve_Freq'};

% Safety check
assert(size(spikeanalysis,2) == numel(colHeaders), 'Header count must match number of columns.');

% Write header in first row
writecell(colHeaders, spikeanalysis_outXlsx, 'FileType','spreadsheet');

% Write the numeric data starting at A2
writematrix(spikeanalysis, spikeanalysis_outXlsx, 'FileType','spreadsheet', 'Range','A2');

