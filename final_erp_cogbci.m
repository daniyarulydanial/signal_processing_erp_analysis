%% USER SETTINGS (edit these paths / flags as needed)
data_root = 'C:/Users/User/Jupyter projects/eeg-analysis/cog_bci';
subjects = {'sub-01', 'sub-03', 'sub-04'};           % subject IDs
sessions = {'ses-01'};           % session IDs
task_filename = 'twoBACK.set';                     % filename to test for each session (case-sensitive)
triggerlist_path = fullfile(data_root,'triggerlist.txt');

OUT_ROOT = fullfile(data_root,'results_paperPipeline'); % outputs go here
PROCESS_SINGLE = false;   % run only the first subject/session (quick test)
PROCESS_ALL    = true;  % set to true to run whole dataset (overrides PROCESS_SINGLE)
MIN_TRIALS_PER_BIN = 30;  % warn/flag if fewer trials

% Preproc params (paper)
DOWNSAMPLE_TO = 250;     % Hz
HP_CUTOFF = 1;           % Hz (high-pass)
CLEANLINE_FREQ = 50;     % Hz
ICLABEL_THRESHOLD = 0.90; % threshold for Eye/Muscle/Heart
EPOCH_MS = [-200 800];   % epoch window (ms)
BASELINE_MS = [-200 0];  % baseline (ms)
AMP_THRESH_UV = 100;     % amplitude rejection ±uV
MWP_WINDOW = 200;        % moving-window width (ms) for pop_artmwppth
MWP_STEP = 50;           % moving-window step (ms)
CHANNEL_SD_CRIT = 2;     % channels with SD > mean+2*sd are bad

% Save options
SAVE_INTERMEDIATE = true;
PLOT_PNGS = true;

%% Setup output folders and logfile
if ~exist(OUT_ROOT,'dir'), mkdir(OUT_ROOT); end
proc_folder = fullfile(OUT_ROOT,'processed'); if ~exist(proc_folder,'dir'), mkdir(proc_folder); end
erp_folder  = fullfile(OUT_ROOT,'erp');       if ~exist(erp_folder,'dir'), mkdir(erp_folder); end
plots_folder= fullfile(OUT_ROOT,'plots');     if ~exist(plots_folder,'dir'), mkdir(plots_folder); end
logfile = fullfile(OUT_ROOT,'pipeline_log.txt');
fidlog = fopen(logfile,'w');
if fidlog==-1, error('Cannot open log file %s', logfile); end
fprintf(fidlog, 'Pipeline started: %s\n', datestr(now));

% --- Prevent figure windows from popping up (batch/headless mode) ---
set(0, 'DefaultFigureVisible', 'off');                        
set(0, 'DefaultFigureCreateFcn', @(fig,~) set(fig, 'Visible', 'off')); 
set(0, 'DefaultFigureToolbar', 'none');
set(0, 'DefaultFigureNumberTitle', 'off');
set(groot,'defaultTextInterpreter','none');
set(groot,'defaultAxesTickLabelInterpreter','none');
set(groot,'defaultLegendInterpreter','none');

%% Start EEGLAB (nogui)
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui');
fprintf(fidlog, 'EEGLAB started: %s\n', datestr(now));

%% Read triggerlist (CSV-like) -> map code->label
fprintf('Parsing triggerlist: %s\n', triggerlist_path);
fprintf(fidlog, 'Parsing triggerlist: %s\n', triggerlist_path);
if ~exist(triggerlist_path,'file')
    error('triggerlist.txt not found at %s', triggerlist_path);
end
T = readtable(triggerlist_path,'Delimiter',',','ReadVariableNames',true,'TextType','string');
codes = cellstr(strtrim(string(T.code)));
labels = cellstr(strtrim(string(T.content)));
labels = regexprep(labels,';',''); % remove stray semicolons
trig_map = containers.Map;
for i=1:numel(codes)
    trig_map(codes{i}) = labels{i};
end
fprintf(fidlog, 'Loaded %d trigger mappings\n', numel(codes));


% ----------------------------
% Helper: robust read of ntrials for bin b
% ----------------------------
function ntr = get_bin_ntrials(ERP, b)
    % Robustly extract *accepted* trial count for bin b from various ERP.ntrials encodings.
    ntr = NaN;
    if ~exist('ERP','var') || isempty(ERP)
        return;
    end
    if isfield(ERP,'ntrials') && ~isempty(ERP.ntrials)
        val = ERP.ntrials;
        try
            % Case 1: struct with field 'accepted' (common ERPLAB)
            if isstruct(val)
                if isfield(val,'accepted')
                    ac = val.accepted;
                    if isnumeric(ac)
                        if numel(ac) >= b, ntr = double(ac(b)); end
                    elseif iscell(ac)
                        try ntr = double(ac{min(numel(ac),b)}); catch, ntr = str2double(char(ac{min(numel(ac),b)})); end
                    end
                    return;
                end
                % Other possible numeric fields (try to pick a numeric vector)
                fn = fieldnames(val);
                for k = 1:numel(fn)
                    tmp = val.(fn{k});
                    if isnumeric(tmp) && numel(tmp) >= b
                        ntr = double(tmp(b)); return;
                    elseif iscell(tmp) && numel(tmp) >= b
                        try ntr = double(tmp{b}); return; catch, ntr = str2double(char(tmp{b})); return; end
                    end
                end
            end

            % Case 2: plain numeric vector
            if isnumeric(val)
                if numel(val) >= b, ntr = double(val(b));
                elseif numel(val) == 1, ntr = double(val(1)); end
                return;
            end

            % Case 3: cell array of values/strings
            if iscell(val)
                if numel(val) >= b
                    v = val{b};
                    if isnumeric(v), ntr = double(v);
                    else ntr = str2double(char(v)); end
                elseif ~isempty(val)
                    v = val{1};
                    if isnumeric(v), ntr = double(v); else ntr = str2double(char(v)); end
                end
                return;
            end

            % Case 4: string/char containing numbers
            if ischar(val) || isstring(val)
                s = char(val);
                toks = regexp(s,'[\d]+','match');
                if numel(toks) >= b, ntr = str2double(toks{b});
                elseif ~isempty(toks), ntr = str2double(toks{1}); end
                return;
            end

        catch
            % fall through to NaN
            ntr = NaN;
            return;
        end
    end

    % Final fallbacks: check other possible fields used by some pipelines
    try
        if isfield(ERP,'ntrials_vec') && numel(ERP.ntrials_vec) >= b
            ntr = double(ERP.ntrials_vec(b)); return;
        end
        if isfield(ERP,'ntrials_per_bin') && numel(ERP.ntrials_per_bin) >= b
            ntr = double(ERP.ntrials_per_bin(b)); return;
        end
    catch
    end
end

% ----------------------------
% Helper: check all bins and produce boolean list of bins meeting threshold
% ----------------------------
function ok = bins_with_min_trials(ERP, min_trials)
    if ~isfield(ERP,'nbin') || isempty(ERP.nbin)
        if isfield(ERP,'bindata'), nbin = size(ERP.bindata,3); else nbin = 0; end
    else
        nbin = ERP.nbin;
    end
    ok = false(1, nbin);
    for bb = 1:nbin
        ntr = get_bin_ntrials(ERP, bb);
        ok(bb) = ~isnan(ntr) && (ntr >= min_trials);
    end
end

%% Helper: write bindesc file for given task (filter triggers that mention TWOBACK)
function bdf_path = write_bindesc_for_task(outdir, trig_map, task_keyword, subj_tag)
    if nargin<4, subj_tag = ''; end
    bdf_path = fullfile(outdir, sprintf('bindesc_%s%s.txt', lower(task_keyword), subj_tag));
    fid = fopen(bdf_path,'w');
    if fid==-1, error('Cannot write bindesc to %s', bdf_path); end
    keys_ = keys(trig_map);
    binidx = 1;
    for k=1:numel(keys_)
        lab = trig_map(keys_{k});
        if contains(lab, upper(task_keyword), 'IgnoreCase', true)
            code = keys_{k};
            fprintf(fid, 'bin %d\n', binidx);
            safe_label = matlab.lang.makeValidName(strrep(lab,' ','_'));
            fprintf(fid, '%s\n', safe_label);
            fprintf(fid, '.{%s}\n\n', code);
            binidx = binidx + 1;
        end
    end
    fclose(fid);
end

%% Processing core: function to process one file (subject, session)
function summary = process_one(subj, ses, data_root, proc_folder, erp_folder, plots_folder, trig_map, triggerlist_path, ...
        DOWNSAMPLE_TO, HP_CUTOFF, CLEANLINE_FREQ, CHANNEL_SD_CRIT, MIN_TRIALS_PER_BIN, IC_LABEL_THRESH, ...
        EPOCH_MS, BASELINE_MS, AMP_THRESH_UV, MWP_WINDOW, MWP_STEP, SAVE_INTERMEDIATE, PLOT_PNGS, fidlog)
    summary = struct();
    try
        subj_eeg_dir = fullfile(data_root, subj, ses, 'eeg');
        listing = dir(fullfile(subj_eeg_dir,'*two*BACK*.set'));
        if isempty(listing)
            listing = dir(fullfile(subj_eeg_dir,'twoBACK.set'));
        end
        if isempty(listing)
            error('twoBACK .set not found for %s/%s in %s', subj, ses, subj_eeg_dir);
        end
        setfile = fullfile(listing(1).folder, listing(1).name);
        fprintf('Processing %s\n', setfile);
        fprintf(fidlog, '[%s] Processing %s\n', datestr(now), setfile);
        
        % Load dataset
        EEG = pop_loadset('filename', listing(1).name, 'filepath', listing(1).folder);
        EEG = eeg_checkset(EEG);
        summary.orig_nbchan = EEG.nbchan;
        summary.srate_orig = EEG.srate;
        
        % --- SAVE ORIGINAL CHANLOCS AND PRECOMPUTE SCALP/PHYSIO LABELS ---
        full_chanlocs = EEG.chanlocs;                     
        full_chan_labels = {full_chanlocs.labels};
        pre_nonEEG_patterns = {'ECG','EOG','EMG','nas','lhj','rhj','EXG'};
        nonEEG_mask_full = false(1,numel(full_chan_labels));
        for p = 1:numel(pre_nonEEG_patterns)
            nonEEG_mask_full = nonEEG_mask_full | contains(full_chan_labels, pre_nonEEG_patterns{p}, 'IgnoreCase', true);
        end
        nonEEG_labels_full = full_chan_labels(nonEEG_mask_full);     
        scalp_chanlocs = full_chanlocs(~nonEEG_mask_full);           
        % --- end precompute ---
        
        % try load chanlocs file in session/chanlocs/get_chanlocs.txt
        chfile = fullfile(data_root, subj, ses, 'chanlocs', 'get_chanlocs.txt');
        if exist(chfile,'file')
            fprintf(fidlog, 'Found chanlocs file: %s\n', chfile);
            fid_ch = fopen(chfile,'r');
            if fid_ch == -1
                fprintf(fidlog, 'Could not open chanlocs file: %s\n', chfile);
                parsed = {};
            else
                parsed = {};
                tline = fgetl(fid_ch);
                while ischar(tline)
                    tline = strtrim(tline);
                    if isempty(tline)
                        tline = fgetl(fid_ch);
                        continue;
                    end
                    parts = regexp(tline, '\s+', 'split');
                    if numel(parts) >= 4
                        lab = parts{1};
                        x = str2double(parts{2});
                        y = str2double(parts{3});
                        z = str2double(parts{4});
                        if ~(isnan(x) || isnan(y) || isnan(z))
                            parsed(end+1,:) = {char(lab), x, y, z}; %#ok<SAGROW>
                        end
                    end
                    tline = fgetl(fid_ch);
                end
                fclose(fid_ch);
            end
            
            if ~isempty(parsed)
                chanlabs = {EEG.chanlocs.labels};
                nParsed = size(parsed,1);
                matched = 0;
                for p=1:nParsed
                    lab = parsed{p,1};
                    idx = find(strcmpi(chanlabs, lab),1);
                    if ~isempty(idx)
                        EEG.chanlocs(idx).X = parsed{p,2};
                        EEG.chanlocs(idx).Y = parsed{p,3};
                        EEG.chanlocs(idx).Z = parsed{p,4};
                        matched = matched + 1;
                    end
                end
                fprintf(fidlog, 'Assigned coordinates for %d channels based on %s\n', matched, chfile);
                summary.chanlocs_assigned = matched;
            else
                fprintf(fidlog, 'Parsed chanlocs file but no coordinates found or file empty.\n');
                summary.chanlocs_assigned = 0;
            end
        end
        
        %% Downsample
        if EEG.srate > DOWNSAMPLE_TO
            EEG = pop_resample(EEG, DOWNSAMPLE_TO);
            EEG = eeg_checkset(EEG);
            fprintf(fidlog, 'Downsampled to %d Hz\n', DOWNSAMPLE_TO);
        else
            fprintf(fidlog, 'No downsampling required (srate=%d)\n', EEG.srate);
        end
        
        %% High-pass filter (1 Hz)
        EEG = pop_eegfiltnew(EEG, 'locutoff', HP_CUTOFF, 'plotfreqz',0);
        EEG = eeg_checkset(EEG);
        fprintf(fidlog, 'High-pass filtered at %0.2f Hz\n', HP_CUTOFF);
        
        %% Cleanline (50 Hz) if available
        if exist('pop_cleanline','file')==2
            try
                EEG = pop_cleanline(EEG,'Bandwidth',2,'ChanCompIndices',1:EEG.nbchan,...
                    'LineFrequencies',CLEANLINE_FREQ,'NormalizeSpectrum','on','PlotFigures','off');
                fprintf(fidlog, 'Applied CleanLine at %d Hz\n', CLEANLINE_FREQ);
            catch ME
                fprintf(fidlog, 'Cleanline failed: %s\n', ME.message);
            end
        else
            fprintf(fidlog, 'cleanline plugin not found, skipping CleanLine\n');
        end
        
        % Save pre-ICA copy if desired
        preica_name = sprintf('%s_%s_preICA.set', subj, ses);
        if SAVE_INTERMEDIATE
            EEG = pop_saveset(EEG, 'filename', preica_name, 'filepath', proc_folder);
            fprintf(fidlog,'Saved pre-ICA set: %s\n', fullfile(proc_folder,preica_name));
        end
        
        % ------------------- Remove non-EEG channels for ICA/bad-channel detection -------------------
        allLabels = {EEG.chanlocs.labels};
        nonEEG_labels = nonEEG_labels_full;
        nonEEG_idx = nonEEG_mask_full;
        
        if ~isempty(nonEEG_labels)
            try
                physioEEG = pop_select(EEG, 'channel', nonEEG_labels);
                physio_name = sprintf('%s_%s_physioChannels.set', subj, ses);
                pop_saveset(physioEEG, 'filename', physio_name, 'filepath', proc_folder);
                fprintf(fidlog, 'Saved physio-only set: %s\n', fullfile(proc_folder, physio_name));
            catch ME
                fprintf(fidlog, 'Warning: could not save physio-only set: %s\n', ME.message);
            end
        
            try
                EEG = pop_select(EEG, 'nochannel', nonEEG_labels);
                EEG = eeg_checkset(EEG);
                fprintf(fidlog, 'Excluded non-EEG channels for ICA/bad-channel detection: %s\n', strjoin(nonEEG_labels,', '));
            catch ME
                fprintf(fidlog, 'Error removing non-EEG channels: %s\n', ME.message);
            end
        else
            fprintf(fidlog, 'No non-EEG channels detected (continuing with full montage)\n');
        end
        
        %% Automatic bad-channel detection: 2-SD rule on channel std
        chan_std = std(double(EEG.data),0,2);
        mean_std = mean(chan_std);
        sd_std = std(chan_std);
        bad_thresh_high = mean_std + CHANNEL_SD_CRIT*sd_std;
        bad_thresh_low = mean_std - CHANNEL_SD_CRIT*sd_std;
        badChanIdx = find(chan_std > bad_thresh_high | chan_std < bad_thresh_low);
        badChanLabels = {EEG.chanlocs(badChanIdx).labels};
        summary.bad_channels = badChanLabels;
        fprintf(fidlog,'Auto-detected %d bad channels by 2-SD rule: %s\n', numel(badChanIdx), strjoin(badChanLabels,', '));
        
        % Interpolate bad channels
        if ~isempty(badChanIdx)
            try
                badLabels = {EEG.chanlocs(badChanIdx).labels};
                EEG = pop_select(EEG, 'nochannel', badLabels);
                EEG = eeg_checkset(EEG);
                EEG = pop_interp(EEG, scalp_chanlocs, 'spherical'); 
                EEG = eeg_checkset(EEG);
                fprintf(fidlog,'Interpolated %d channels and restored to scalp montage\n', numel(badChanIdx));
            catch ME
                fprintf(fidlog,'Interpolation failed: %s\n', ME.message);
            end
        else
            fprintf(fidlog,'No channels interpolated.\n');
        end
        summary.num_interpolated = numel(badChanIdx);
        summary.interpolated_labels = badChanLabels;
        
        %% Full-rank average reference
        EEG = pop_reref(EEG, []);
        EEG = eeg_checkset(EEG);
        fprintf(fidlog,'Re-referenced to average (full-rank)\n');
        
        %% Run ICA (runica extended)
        fprintf(fidlog,'Running ICA (runica extended)... this may take a while\n');
        try
            EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1, 'interrupt','off');
            EEG = eeg_checkset(EEG);
            fprintf(fidlog,'ICA finished: %d components\n', size(EEG.icaweights,1));
        catch ME
            error('ICA failed: %s', ME.message);
        end
        
        postica_name = sprintf('%s_%s_postICA.set', subj, ses);
        if SAVE_INTERMEDIATE
            EEG = pop_saveset(EEG, 'filename', postica_name, 'filepath', proc_folder);
            fprintf(fidlog,'Saved post-ICA set: %s\n', fullfile(proc_folder,postica_name));
        end
        
        %% ICLabel classification and automatic component removal (robust save + removal)
        removedICs = [];
        if exist('pop_iclabel','file')==2
            EEG = pop_iclabel(EEG,'default');
            probs = EEG.etc.ic_classification.ICLabel.classifications;
            nc = size(probs,1);
            ncols = size(probs,2);
            if ncols == 6
                colnames = {'Brain','Muscle','Eye','Heart','LineNoise','Other'};
            elseif ncols == 7
                colnames = {'Brain','Muscle','Eye','Heart','LineNoise','Other','Unknown'};
            else
                colnames = arrayfun(@(x) sprintf('p%d',x), 1:ncols, 'UniformOutput', false);
            end
        
            try
                ICidx = (1:nc);
                T_all = array2table([ICidx, probs], 'VariableNames', [{'IC'}, colnames]);
                writetable(T_all, fullfile(proc_folder, sprintf('%s_%s_ICLabel_all.csv', subj, ses)));
                fprintf(fidlog, 'Saved ICLabel table (all ICs): %s\n', fullfile(proc_folder, sprintf('%s_%s_ICLabel_all.csv', subj, ses)));
            catch ME
                fprintf(fidlog, 'Failed to save full ICLabel CSV: %s\n', ME.message);
            end
        
            eye_idx   = find(probs(:, strcmpi(colnames, 'Eye'))   >= IC_LABEL_THRESH);
            mus_idx   = find(probs(:, strcmpi(colnames, 'Muscle'))>= IC_LABEL_THRESH);
            heart_idx = [];
            if any(strcmpi(colnames,'Heart'))
                heart_idx = find(probs(:, strcmpi(colnames,'Heart')) >= IC_LABEL_THRESH);
            end
            removedICs = unique([eye_idx; mus_idx; heart_idx]);
        
            try
                if ~isempty(removedICs)
                    removedMat = [removedICs, probs(removedICs,:)];
                    T_removed = array2table(removedMat, 'VariableNames', [{'IC'}, colnames]);
                    outCSV = fullfile(proc_folder, sprintf('%s_%s_removedICs.csv', subj, ses));
                    writetable(T_removed, outCSV);
                    fprintf(fidlog, 'Saved removed ICs CSV: %s (rows=%d)\n', outCSV, size(T_removed,1));
                else
                    fprintf(fidlog, 'No ICs exceeded ICLabel threshold %.2f\n', IC_LABEL_THRESH);
                end
            catch ME
                fprintf(fidlog, 'Failed to save removed ICs CSV: %s\n', ME.message);
            end

            try
                if PLOT_PNGS && ~isempty(removedICs)
                    figdir = fullfile(plots_folder, subj, ses, 'ICA'); 
                    if ~exist(figdir,'dir')
                        mkdir(figdir); 
                    end

                    save_ic_topoplot(EEG, removedICs, figdir, subj, ses);
                    save_ic_spectra(EEG, removedICs, figdir, subj, ses);
                    fprintf(fidlog, 'Saved IC figures for %d removed ICs\n', numel(removedICs));
                end
            catch ME
                fprintf(fidlog, 'Failed to save IC figures: %s\n', ME.message); 
            end
        
            if ~isempty(removedICs)
                EEG = pop_subcomp(EEG, removedICs, 0);
                EEG = eeg_checkset(EEG);
                fprintf(fidlog, 'Removed %d ICs: %s\n', numel(removedICs), mat2str(removedICs));
            end
        
        else
            fprintf(fidlog, 'ICLabel not found - skipping automatic IC removal (please inspect manually)\n');
        end

        cleaned_name = sprintf('%s_%s_cleaned.set', subj, ses);
        if SAVE_INTERMEDIATE
            EEG = pop_saveset(EEG, 'filename', cleaned_name, 'filepath', proc_folder);
            fprintf(fidlog,'Saved ICA-cleaned dataset: %s\n', fullfile(proc_folder,cleaned_name));
        end
        
        try
            drawnow limitrate; 
            close all hidden; 
            pause(0.05); 
        catch
        end

        %% Create ERPLAB eventlist (required for binlister)
        try
            EEG = pop_creabasiceventlist(EEG, 'AlphanumericCleaning','on','Newboundary',{-99}, 'Stringboundary',{'boundary'});
            fprintf(fidlog,'Created basic ERPLAB eventlist\n');
        catch
            fprintf(fidlog,'pop_creabasiceventlist failed - continuing but ERPLAB binlister may fail\n');
        end
        
        %% Create bindesc file for TWOBACK triggers only
        bindesc_file = write_bindesc_for_task(proc_folder, trig_map, 'TWOBACK', sprintf('_%s_%s', subj, ses));
        fprintf(fidlog,'Wrote bindesc file: %s\n', bindesc_file);
        
        %% Run Binlister (assign bins)
        try
            export_el_file = fullfile(proc_folder, sprintf('%s_%s_eventlist_binned.txt', subj, ses));
            EEG = pop_binlister(EEG, 'BDF', bindesc_file, 'IndexEL',1, 'ExportEL', export_el_file, 'SendEL2','EEG','UpdateEEG','on', 'Warning','on');
            fprintf(fidlog,'Ran pop_binlister; exported eventlist to %s\n', export_el_file);
        catch ME
            fprintf(fidlog,'pop_binlister failed: %s\n', ME.message);
        end
        
        %% Epoch bins (ERPLAB pop_epochbin)
        try
            EEG = pop_epochbin(EEG, EPOCH_MS, BASELINE_MS);
            EEG = eeg_checkset(EEG);
            fprintf(fidlog,'Epoching done [%d %d] ms with baseline [%d %d] ms\n', EPOCH_MS(1), EPOCH_MS(2), BASELINE_MS(1), BASELINE_MS(2));
        catch ME
            error('Epoching failed: %s', ME.message);
        end
        
        % Save epoched set
        epoched_name = sprintf('%s_%s_epoched.set', subj, ses);
        if SAVE_INTERMEDIATE
            EEG = pop_saveset(EEG, 'filename', epoched_name, 'filepath', proc_folder);
            fprintf(fidlog,'Saved epoched set: %s\n', fullfile(proc_folder, epoched_name));
        end
        
        try
            drawnow limitrate; 
            close all hidden; 
            pause(0.05); 
        catch
        end

        %% Artifact detection on epochs
        try
            EEG = pop_artextval(EEG, 'Channel', 1:EEG.nbchan, 'Flag', 1, 'Threshold', [-AMP_THRESH_UV AMP_THRESH_UV], 'Twindow', EPOCH_MS);
            fprintf(fidlog, 'Applied amplitude threshold ±%d µV\n', AMP_THRESH_UV);
        catch ME
            fprintf(fidlog, 'pop_artextval failed: %s\n', ME.message);
        end
        try
            EEG = pop_artmwppth(EEG, 'Channel', 1:EEG.nbchan, 'Flag', [1 4], 'Review', 'off', 'Threshold', AMP_THRESH_UV, 'Twindow', EPOCH_MS, 'Windowsize', MWP_WINDOW, 'Windowstep', MWP_STEP);
            fprintf(fidlog, 'Applied moving-window (win=%d ms, step=%d ms) peak-to-peak threshold %d µV\n', MWP_WINDOW, MWP_STEP, AMP_THRESH_UV);
        catch ME
            fprintf(fidlog, 'pop_artmwppth failed: %s\n', ME.message);
        end
        
        try
            drawnow limitrate; 
            close all hidden; 
            pause(0.05); 
        catch
        end

        %% Epoch counts summary (existing code unchanged)
        try
            nEpochs = EEG.trials;
            epoch_bin = cell(nEpochs,1);
            epoch_rej = false(nEpochs,1);
            for e = 1:nEpochs
                if isfield(EEG.epoch(e), 'eventtype')
                    evt = EEG.epoch(e).eventtype;
                    if iscell(evt), epoch_bin{e} = evt{1}; else epoch_bin{e} = evt; end
                else
                    if ~isempty(EEG.epoch(e).event)
                        if isstruct(EEG.epoch(e).event)
                            epoch_bin{e} = EEG.epoch(e).event(1).type;
                        elseif iscell(EEG.epoch(e).event)
                            epoch_bin{e} = EEG.epoch(e).event{1}.type;
                        else
                            epoch_bin{e} = '';
                        end
                    else
                        epoch_bin{e} = '';
                    end
                end
                if isfield(EEG.reject, 'rejglobal') && ~isempty(EEG.reject.rejglobal)
                    try
                        epoch_rej(e) = any(EEG.reject.rejglobal(:,e));
                    catch
                        epoch_rej(e) = false;
                    end
                else
                    if isfield(EEG.reject,'rejmanual')
                        try, epoch_rej(e) = any(EEG.reject.rejmanual(:,e)); catch, end
                    end
                end
            end
            uniq_bins = unique(epoch_bin);
            bin_summary = table;
            bin_labels = {};
            accepted_counts = [];
            total_counts = [];
            for b=1:numel(uniq_bins)
                thisbin = uniq_bins{b};
                idx = find(cellfun(@(x) isequal(x,thisbin), epoch_bin));
                total_counts(b) = numel(idx);
                accepted_counts(b) = sum(~epoch_rej(idx));
                bin_labels{b} = thisbin;
            end
            bin_summary.Bin = bin_labels;
            bin_summary.Total = total_counts;
            bin_summary.Accepted = accepted_counts;
            writetable(bin_summary, fullfile(proc_folder, sprintf('%s_%s_epochCounts.csv', subj, ses)));
            fprintf(fidlog,'Saved epoch counts per bin CSV\n');
            summary.bin_summary = bin_summary;

            try
                if PLOT_PNGS
                    figdir = fullfile(plots_folder, subj, ses); 
                    if ~exist(figdir,'dir'), mkdir(figdir); end
                    save_epoch_counts(bin_summary, figdir, subj, ses);
                    fprintf(fidlog,'Saved epoch counts figure\n');
                    summary.bin_summary = bin_summary;
                end
            catch ME
                fprintf(fidlog,'Failed to save epoch counts figure: %s\n', ME.message); 
            end

        catch ME
            fprintf(fidlog,'Failed to produce epoch counts: %s\n', ME.message);
        end
        
        try
            drawnow limitrate; 
            close all hidden; 
            pause(0.05); 
        catch
        end

        %% Averaging -> ERP (ERPLAB)
        % === PATCHED: split averaging call from post-averaging bookkeeping for accurate error reporting ===
        try
            ERP = pop_averager(EEG, 'Criterion', 'good', 'DSindex', 1, 'ExcludeBoundary', 'on', 'SEM', 'on');
            ERP.erpname = sprintf('%s_%s_ERP', subj, ses);
        catch ME
            fprintf(fidlog,'pop_averager failed: %s\n', ME.message);
            % skip post-averaging steps for this subject/session
            summary.status = 'averaging_failed';
        end

        % Post-averaging bookkeeping (separate try/catch so errors here are not attributed to pop_averager)
        try
            % --- determine per-bin trial counts and which bins meet MIN_TRIALS_PER_BIN ---
            % robustly determine number of bins
            if isfield(ERP,'nbin') && ~isempty(ERP.nbin)
                nbin = double(ERP.nbin);
            elseif isfield(ERP,'bindata') && ~isempty(ERP.bindata)
                nbin = size(ERP.bindata, 3);
            else
                nbin = 0;
            end
            nbin = max(0, nbin);
            
            ntrials_vec = nan(1, max(1,nbin));
            for bb = 1:nbin
                ntrials_vec(bb) = get_bin_ntrials(ERP, bb);
            end
            
            % Define validity using explicit numeric check & threshold
            valid_bins = (~isnan(ntrials_vec)) & (ntrials_vec >= MIN_TRIALS_PER_BIN);
            
            % Save per-bin trial info (CSV)
            try
                Tper = table((1:nbin), ntrials_vec', valid_bins', 'VariableNames', {'bin','ntrials','meets_min'});
                writetable(Tper, fullfile(proc_folder, sprintf('%s_%s_perbin_trials.csv', subj, ses)));
                fprintf(fidlog,'Saved per-bin trial counts CSV (%d bins). Counts: %s\n', nbin, mat2str(ntrials_vec));
                fprintf('[%s %s] Bins meeting MIN_TRIALS=%d : %s\n', subj, ses, MIN_TRIALS_PER_BIN, mat2str(find(valid_bins)));
            catch ME
                fprintf(fidlog,'Failed to write per-bin CSV: %s\n', ME.message);
            end

            % === PATCHED: ensure ERP.ntrials.accepted exists for downstream scripts ===
            try
                if ~isfield(ERP,'ntrials') || isempty(ERP.ntrials)
                    ERP.ntrials = struct();
                end
                ERP.ntrials.accepted = ntrials_vec;
            catch
                fprintf(fidlog,'Warning: Could not assign ERP.ntrials.accepted\n');
            end

            % Save ERP (always) and print human-friendly bins passing message
            try
                pop_savemyerp(ERP, 'erpname', ERP.erpname, 'filename', [ERP.erpname '.erp'], 'filepath', erp_folder, 'Warning', 'off');
            catch ME
                fprintf(fidlog,'Failed to save ERP file: %s\n', ME.message);
            end

            if all(valid_bins)
                fprintf(fidlog,'Saved ERP: %s\n', fullfile(erp_folder,[ERP.erpname '.erp']));
            else
                fprintf(fidlog,'Saved ERP (but some bins below %d trials): %s\n', MIN_TRIALS_PER_BIN, fullfile(erp_folder,[ERP.erpname '.erp']));
                bins_passing = find(valid_bins);
                if isempty(bins_passing)
                    fprintf(fidlog,'No bins pass threshold %d.\n', MIN_TRIALS_PER_BIN);
                else
                    fprintf(fidlog,'Bins passing threshold: %s\n', mat2str(bins_passing));
                end
                ERP.user_min_trials_ok = valid_bins;
            end

            % Save ERP figures (traces + topomaps) - keep this in its own try/catch
            try
                if PLOT_PNGS
                    figdir = fullfile(plots_folder, subj, ses, 'ERP'); 
                    if ~exist(figdir,'dir'), mkdir(figdir); end
                    % standard channels
                    save_erp_traces(ERP, {'Cz','Pz','Fz'}, figdir, subj, ses);
                    % topomap windows (example P300 window)
                    save_erp_topomap(ERP, [250 450], figdir, subj, ses, 'P300');
                    fprintf(fidlog,'Saved ERP figures\n');
                end
            catch ME 
                fprintf(fidlog,'Failed to save ERP figures: %s\n', ME.message); 
            end

        catch ME
            fprintf(fidlog,'ERP post-processing failed: %s\n', ME.message);
        end

        try 
            drawnow limitrate; 
            close all hidden; 
            pause(0.05); 
        catch 
        end
               
        % Return summary struct
        summary.status = 'ok';
        summary.setfile = setfile;
        
        if exist('ERP','var') && isstruct(ERP) && isfield(ERP,'erpname')
            summary.ERPname = ERP.erpname;
        end

    catch ME
        summary.status = 'error';
        summary.message = ME.message;
        fprintf(fidlog,'ERROR processing %s %s: %s\n', subj, ses, ME.message);
    end
end

%% MAIN: run single test or full processing
if PROCESS_ALL
    PROCESS_SINGLE = false;
end

if PROCESS_SINGLE
    s0 = subjects{1};
    se0 = sessions{1};
    fprintf('=== Running single test on %s/%s ===\n', s0, se0);
    fprintf(fidlog,'=== Running single test on %s/%s ===\n', s0, se0);
    summary1 = process_one(s0, se0, data_root, proc_folder, erp_folder, plots_folder, trig_map, triggerlist_path, ...
        DOWNSAMPLE_TO, HP_CUTOFF, CLEANLINE_FREQ, CHANNEL_SD_CRIT, MIN_TRIALS_PER_BIN, ICLABEL_THRESHOLD, ...
        EPOCH_MS, BASELINE_MS, AMP_THRESH_UV, MWP_WINDOW, MWP_STEP, SAVE_INTERMEDIATE, PLOT_PNGS, fidlog);
    disp(summary1);
end

if PROCESS_ALL
    ALLERP = buildERPstruct([]); CURRENTERP = 0;
    ERP_list = {};
    for si = 1:numel(subjects)
        for se = 1:numel(sessions)
            subj = subjects{si};
            ses = sessions{se};
            fprintf('Processing %s/%s ...\n', subj, ses);
            ssum = process_one(subj, ses, data_root, proc_folder, erp_folder, plots_folder, trig_map, triggerlist_path, ...
                DOWNSAMPLE_TO, HP_CUTOFF, CLEANLINE_FREQ, CHANNEL_SD_CRIT, MIN_TRIALS_PER_BIN, ICLABEL_THRESHOLD, ...
                EPOCH_MS, BASELINE_MS, AMP_THRESH_UV, MWP_WINDOW, MWP_STEP, SAVE_INTERMEDIATE, PLOT_PNGS, fidlog);
            if isfield(ssum,'ERPname')
                try
                    ERP = pop_loaderp('filename', [ssum.ERPname '.erp'], 'filepath', erp_folder);
                    CURRENTERP = CURRENTERP + 1;
                    ALLERP(CURRENTERP) = ERP;
                    ERP_list{end+1} = ERP;
                catch
                    fprintf(fidlog, 'Failed to load ERP for %s/%s\n', subj, ses);
                end
            end
        end
    end
    if CURRENTERP > 0
        ERP_G = pop_gaverager(ALLERP, 'Criterion', 100, 'ERPindex', [1:CURRENTERP]);
        ERP_G.erpname = 'grand_avg_all';
        pop_savemyerp(ERP_G, 'erpname', ERP_G.erpname, 'filename', [ERP_G.erpname '.erp'], 'filepath', erp_folder, 'warning', 'off');
        fprintf(fidlog,'Saved grand-average ERP: %s\n', fullfile(erp_folder, [ERP_G.erpname '.erp']));
    end
end

fprintf(fidlog, 'Pipeline finished at %s\n', datestr(now));
fclose(fidlog);
fprintf('Pipeline finished. Log saved to %s\n', logfile);

%% ----------------------- Local plotting helper functions -----------------------
function save_ic_topoplot(EEG, comps, outdir, subj, ses)
    for c = comps(:)
        try
            fig = figure('Visible','off','Units','pixels','Position',[100 100 600 480]);
            if isfield(EEG,'icawinv') && ~isempty(EEG.icawinv)
                topoplot(EEG.icawinv(:,c), EEG.chanlocs, 'electrodes','on');
            else
                text(0.3,0.5,sprintf('No icawinv for comp %d',c));
            end
            title(sprintf('%s %s - IC %d', subj, ses, c),'Interpreter','none');
            fname = fullfile(outdir, sprintf('%s_%s_IC%d_topo.png', subj, ses, c));
            try 
                exportgraphics(fig, fname, 'Resolution',300); 
            catch 
                saveas(fig,fname); 
            end
            close(fig);
        catch 
            try 
                close(fig); 
            catch 
            end 
        end
    end
end

function save_ic_spectra(EEG, comps, outdir, subj, ses)
    for c = comps(:)
        try
            fig = figure('Visible','off','Units','pixels','Position',[100 100 700 450]);
            if isfield(EEG,'icaact') && ~isempty(EEG.icaact)
                icact = EEG.icaact(c,:);
            else
                if isfield(EEG,'icawinv') && ~isempty(EEG.icawinv)
                    icact = EEG.icawinv(:,c) * double(EEG.data);
                else
                    icact = [];
                end
            end
            if ~isempty(icact)
                try
                    spectopo(double(icact), 0, EEG.srate, 'plot','off');
                    title(sprintf('%s %s - IC %d spectrum', subj, ses, c));
                catch
                    [Pxx, F] = pwelch(double(icact), [], [], [], EEG.srate);
                    plot(F,10*log10(Pxx)); xlim([0 50]); xlabel('Hz'); ylabel('dB');
                end
            else
                text(0.3,0.5,'No IC activations available');
            end
            fname = fullfile(outdir, sprintf('%s_%s_IC%d_spec.png', subj, ses, c));
            try 
                exportgraphics(fig, fname, 'Resolution',300); catch, saveas(fig,fname); 
            end
            close(fig);
        catch 
            try 
                close(fig); 
            catch 
            end 
        end
    end
end

function save_epoch_counts(bin_summary, outdir, subj, ses)
    fig = figure('Visible','off','Units','pixels','Position',[100 100 800 300]);
    bar(bin_summary.Accepted);
    set(gca,'XTickLabel', bin_summary.Bin, 'XTickLabelRotation',45);
    ylabel('Accepted trials'); title(sprintf('%s %s - epoch counts', subj, ses));
    fname = fullfile(outdir, sprintf('%s_%s_epochCounts.png', subj, ses));
    try 
        exportgraphics(fig,fname,'Resolution',300); 
    catch 
        saveas(fig,fname); 
    end
    close(fig);
end

% === PATCHED: robust save_erp_traces that handles various bindescr formats and avoids dot-indexing errors ===
function save_erp_traces(ERP, chanLabels, outdir, subj, ses)
    if ~isfield(ERP,'times')
        times = 1;
    else
        times = ERP.times; % ms
    end

    % normalized/robust bindescr -> cellstr
    bindescr_list = {};
    try
        if isfield(ERP,'bindescr') && ~isempty(ERP.bindescr)
            if iscell(ERP.bindescr)
                for ii=1:numel(ERP.bindescr)
                    try
                        bindescr_list{ii} = char(ERP.bindescr{ii});
                    catch
                        try bindescr_list{ii} = string(ERP.bindescr{ii}); catch, bindescr_list{ii} = sprintf('bin%d',ii); end
                    end
                end
            elseif ischar(ERP.bindescr) || isstring(ERP.bindescr)
                bindescr_list = {char(ERP.bindescr)};
            elseif isstruct(ERP.bindescr)
                try
                    bindescr_list = arrayfun(@(s) char(s.label), ERP.bindescr, 'UniformOutput', false);
                catch
                    bindescr_list = {};
                end
            else
                bindescr_list = {};
            end
        end
    catch
        bindescr_list = {};
    end

    % fallback default labels if empty
    if isempty(bindescr_list)
        try
            nbin = max(1, ERP.nbin);
        catch
            if isfield(ERP,'bindata'), nbin = size(ERP.bindata,3); else nbin = 1; end
        end
        bindescr_list = arrayfun(@(k) sprintf('bin%d',k), 1:nbin, 'UniformOutput', false);
    end

    for ch = 1:numel(chanLabels)
        chname = chanLabels{ch}; 
        cidx = find(strcmpi({ERP.chanlocs.labels}, chname), 1); 
        if isempty(cidx), continue; end
        fig = figure('Visible','off','Units','pixels','Position',[100 100 900 500]); hold on;
        legend_labels = {};
        for b = 1:max(1,ERP.nbin)
            try
                meanv = [];
                if isnumeric(ERP.bindata)
                    if size(ERP.bindata,1) >= cidx && size(ERP.bindata,3) >= b
                        meanv = squeeze(ERP.bindata(cidx,:,b));
                    end
                elseif iscell(ERP.bindata)
                    if cidx <= size(ERP.bindata,1) && b <= size(ERP.bindata,2)
                        tmp = ERP.bindata{cidx,b};
                        meanv = tmp(:);
                    end
                end
                if isempty(meanv) || all(isnan(meanv)) || all(meanv==0), continue; end
                if isfield(ERP,'sem') && ~isempty(ERP.sem)
                    try
                        semv = squeeze(ERP.sem(cidx,:,b));
                        x = [times, fliplr(times)]; y = [meanv+semv, fliplr(meanv-semv)];
                        hfill = fill(x, y, 'k'); set(hfill,'FaceAlpha',0.15,'EdgeColor','none');
                    catch
                    end
                end
                plot(times, meanv, 'LineWidth', 1.2);
                % legend label: use bindescr_list if available
                if b <= numel(bindescr_list), legend_labels{end+1} = bindescr_list{b}; else legend_labels{end+1} = sprintf('bin%d',b); end
            catch
            end
        end
        xlabel('Time (ms)'); ylabel('\muV');
        title(sprintf('%s %s - %s', subj, ses, chname), 'Interpreter','none');
        % safe legend
        try
            if ~isempty(legend_labels)
                legend(legend_labels, 'Interpreter','none', 'Location','bestoutside');
            end
        catch
            % skip legend if it errors
        end
        fname = fullfile(outdir, sprintf('%s_%s_ERP_%s.png', subj, ses, chname));
        try 
            exportgraphics(fig,fname,'Resolution',300); 
        catch 
            saveas(fig,fname); 
        end
        close(fig);
    end
end

function save_erp_topomap(ERP, timewin, outdir, subj, ses, winname)
    [~, idx1] = min(abs(ERP.times - timewin(1))); [~, idx2] = min(abs(ERP.times - timewin(2)));
    avgmap = mean(ERP.bindata(:, idx1:idx2, :), 2);
    for b = 1:ERP.nbin
        fig = figure('Visible','off','Units','pixels','Position',[100 100 600 480]);
        topo = squeeze(avgmap(:,1,b));
        try
            topoplot(topo, ERP.chanlocs, 'maplimits','absmax','electrodes','on');
        catch
            imagesc(reshape(topo,[],1)); title('Topo placeholder');
        end
        title(sprintf('%s %s - bin %d - %s (%d-%d ms)', subj, ses, b, winname, timewin(1), timewin(2)),'Interpreter','none');
        fname = fullfile(outdir, sprintf('%s_%s_ERP_topo_bin%d_%s.png', subj, ses, b, winname));
        try 
            exportgraphics(fig, fname, 'Resolution',300); 
        catch 
            saveas(fig,fname); 
        end
        close(fig);
    end
end

% End of file