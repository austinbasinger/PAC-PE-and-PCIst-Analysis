%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    SCRIPT 1: PAC_Compute.m                               %%
%%    - Phase frequency: 10 Hz (alpha) at ALL levels        %%
%%    - Amplitude bands: delta, theta, alpha, beta,         %%
%%      low/mid/high gamma                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1) Setup
clear all; 
close all;
clc;

directory = 'C:\Users\Austin\Data\';
PREPROC_Folder = 'EEG_Data\PREPROC\';
PAC_Folder = 'PAC_MultiBand_All\';

if ~exist([directory, PAC_Folder], 'dir')
    mkdir([directory, PAC_Folder]);
end

subjects = [{'P01'},{'P02'},{'P03'},{'P04'},{'P05'},{'P06'},{'P07'},{'P08'},...
    {'P09'},{'P10'},{'P11'},{'P12'},{'P13'},{'P14'},{'P15'},{'P16'},{'P17'}];

raw_labels = [{'P04'},{'P06'},{'P08'},{'P10'},{'P12'},{'P13'},{'P14'},{'P15'},...
    {'P17'},{'P20'},{'P21'},{'P22'},{'P24'},{'P25'},{'P26'},{'P27'},{'P28'}];

%% 2) Parameters
word_duration = 0.32;
words_per_sentence = 4;
sentence_duration = word_duration * words_per_sentence;
num_sentences_per_trial = 12;
trial_start_offset = 1.0;
TARGET_FS = 256;

chans2exclude = {'A1', 'A2', 'COMNT', 'SCALE'};

% Phase frequency: alpha (10 Hz) at ALL levels
PHASE_FREQ = 10;

% Amplitude bands: all 7
amp_bands = struct();
amp_bands.delta = [1 4];
amp_bands.theta = [4 8];
amp_bands.alpha = [8 13];
amp_bands.beta = [13 30];
amp_bands.low_gamma = [30 45];
amp_bands.mid_gamma = [45 60];
amp_bands.high_gamma = [65 80];
band_names = fieldnames(amp_bands);
n_bands = length(band_names);

n_bins = 18;

% Epoch durations (match word boundaries of stimulus exactly)
epoch_durations = struct();
epoch_durations.words = word_duration;                % 320 ms
epoch_durations.subject_phrases = word_duration * 2;  % 640 ms
epoch_durations.sentences = sentence_duration;        % 1280 ms

fprintf('Cycles of 10 Hz phase per epoch:\n');
fprintf('  Words (320 ms):      %.1f cycles\n', epoch_durations.words * PHASE_FREQ);
fprintf('  Phrases (640 ms):    %.1f cycles\n', epoch_durations.subject_phrases * PHASE_FREQ);
fprintf('  Sentences (1280 ms): %.1f cycles\n', epoch_durations.sentences * PHASE_FREQ);
fprintf('\n');

%% 3) Pre-compute filter coefficients
fprintf('Pre-computing filter coefficients...\n');

% Phase filter at 10 Hz
[b_phase, a_phase] = butter(4, PHASE_FREQ/(TARGET_FS/2), 'low');

% Amplitude filters
b_amp = struct();
a_amp = struct();
for b = 1:n_bands
    [b_amp.(band_names{b}), a_amp.(band_names{b})] = ...
        butter(4, amp_bands.(band_names{b})/(TARGET_FS/2), 'bandpass');
end

fprintf('Filter coefficients ready.\n\n');

%% 4) Subject loop
for s = 1:length(subjects)
    
    fprintf('\n=== Processing Subject %s (raw: %s) (%d/%d) ===\n', ...
        subjects{s}, raw_labels{s}, s, length(subjects));
    
    try
        data_file = [directory, PREPROC_Folder, 'DataCleanICA_AR_', char(subjects(s)), '.mat'];
        if ~exist(data_file, 'file')
            warning('File not found for %s. Skipping.', subjects{s});
            continue;
        end
        
        load(data_file);
        
        Fs_original = data.fsample;
        n_trials = length(data.trial);
        fprintf('Data loaded: %d trials, %d Hz\n', n_trials, Fs_original);
        
        if Fs_original ~= TARGET_FS
            fprintf('  Resampling from %d Hz to %d Hz...\n', Fs_original, TARGET_FS);
            for tr = 1:n_trials
                trial_data = data.trial{tr};
                data.trial{tr} = resample(trial_data', TARGET_FS, Fs_original)';
            end
            data.fsample = TARGET_FS;
        end
        Fs = data.fsample;
        fprintf('  Sampling rate: %d Hz\n', Fs);
        
        chan_idx = find(~ismember(data.label, chans2exclude));
        n_chans = length(chan_idx);
        fprintf('  Using %d EEG channels\n', n_chans);
        
        %% Initialize storage
        PAC_all = struct();
        for b = 1:n_bands
            bn = band_names{b};
            PAC_all.([bn '_words']) = [];
            PAC_all.([bn '_subject_phrases']) = [];
            PAC_all.([bn '_sentences']) = [];
        end
        
        %% WORD-LEVEL PAC
        fprintf('  Computing Word-level PAC...\n');
        
        samples_per_word = round(epoch_durations.words * Fs);
        n_word_epochs = n_trials * num_sentences_per_trial * words_per_sentence;
        word_pac_values = zeros(n_bands, n_chans, n_word_epochs);
        epoch_count = 0;
        
        tic;
        for tr = 1:n_trials
            trial_data = data.trial{tr};
            n_samples = size(trial_data, 2);
            
            for sent = 1:num_sentences_per_trial
                for w = 1:words_per_sentence
                    epoch_count = epoch_count + 1;
                    
                    word_idx = (sent-1) * words_per_sentence + w;
                    word_onset = trial_start_offset + (word_idx - 1) * word_duration;
                    
                    onset_sample = round(word_onset * Fs) + 1;
                    end_sample = onset_sample + samples_per_word - 1;
                    
                    if end_sample > n_samples, end_sample = n_samples; end
                    if onset_sample < 1, onset_sample = 1; end
                    
                    epoch_data = trial_data(chan_idx, onset_sample:end_sample);
                    epoch_data = epoch_data - mean(epoch_data, 2);
                    
                    pac_matrix = compute_pac_all_bands(epoch_data, ...
                        b_phase, a_phase, b_amp, a_amp, band_names, n_bins);
                    word_pac_values(:, :, epoch_count) = pac_matrix;
                end
            end
        end
        word_time = toc;
        fprintf('    Word PAC done in %.1f seconds (%d epochs)\n', word_time, epoch_count);
        
        for b = 1:n_bands
            band_pac = squeeze(mean(word_pac_values(b, :, :), 2));
            PAC_all.([band_names{b} '_words']) = mean(band_pac);
        end
        
        %% SUBJECT PHRASE-LEVEL PAC
        fprintf('  Computing Subject Phrase-level PAC...\n');
        
        samples_per_phrase = round(epoch_durations.subject_phrases * Fs);
        n_phrase_epochs = n_trials * num_sentences_per_trial;
        phrase_pac_values = zeros(n_bands, n_chans, n_phrase_epochs);
        epoch_count = 0;
        
        tic;
        for tr = 1:n_trials
            trial_data = data.trial{tr};
            n_samples = size(trial_data, 2);
            
            for sent = 1:num_sentences_per_trial
                epoch_count = epoch_count + 1;
                
                phrase_onset = trial_start_offset + (sent-1) * sentence_duration;
                
                onset_sample = round(phrase_onset * Fs) + 1;
                end_sample = onset_sample + samples_per_phrase - 1;
                
                if end_sample > n_samples, end_sample = n_samples; end
                if onset_sample < 1, onset_sample = 1; end
                
                epoch_data = trial_data(chan_idx, onset_sample:end_sample);
                epoch_data = epoch_data - mean(epoch_data, 2);
                
                pac_matrix = compute_pac_all_bands(epoch_data, ...
                    b_phase, a_phase, b_amp, a_amp, band_names, n_bins);
                phrase_pac_values(:, :, epoch_count) = pac_matrix;
            end
        end
        phrase_time = toc;
        fprintf('    Phrase PAC done in %.1f seconds (%d epochs)\n', phrase_time, epoch_count);
        
        for b = 1:n_bands
            band_pac = squeeze(mean(phrase_pac_values(b, :, :), 2));
            PAC_all.([band_names{b} '_subject_phrases']) = mean(band_pac);
        end
        
        %% SENTENCE-LEVEL PAC
        fprintf('  Computing Sentence-level PAC...\n');
        
        samples_per_sentence = round(epoch_durations.sentences * Fs);
        n_sentence_epochs = n_trials * num_sentences_per_trial;
        sentence_pac_values = zeros(n_bands, n_chans, n_sentence_epochs);
        epoch_count = 0;
        
        tic;
        for tr = 1:n_trials
            trial_data = data.trial{tr};
            n_samples = size(trial_data, 2);
            
            for sent = 1:num_sentences_per_trial
                epoch_count = epoch_count + 1;
                
                sentence_onset = trial_start_offset + (sent-1) * sentence_duration;
                
                onset_sample = round(sentence_onset * Fs) + 1;
                end_sample = onset_sample + samples_per_sentence - 1;
                
                if end_sample > n_samples, end_sample = n_samples; end
                if onset_sample < 1, onset_sample = 1; end
                
                epoch_data = trial_data(chan_idx, onset_sample:end_sample);
                epoch_data = epoch_data - mean(epoch_data, 2);
                
                pac_matrix = compute_pac_all_bands(epoch_data, ...
                    b_phase, a_phase, b_amp, a_amp, band_names, n_bins);
                sentence_pac_values(:, :, epoch_count) = pac_matrix;
            end
        end
        sent_time = toc;
        fprintf('    Sentence PAC done in %.1f seconds (%d epochs)\n', sent_time, epoch_count);
        
        for b = 1:n_bands
            band_pac = squeeze(mean(sentence_pac_values(b, :, :), 2));
            PAC_all.([band_names{b} '_sentences']) = mean(band_pac);
        end
        
        %% Save
        fprintf('  Saving...\n');
        save([directory, PAC_Folder, 'PAC_MultiBand_All_', char(subjects(s)), '.mat'], ...
            'PAC_all', 'Fs', 'PHASE_FREQ', 'amp_bands', 'epoch_durations', '-v7.3');
        
        fprintf('  TOTAL: %.1f seconds for %s\n', word_time + phrase_time + sent_time, subjects{s});
        
        fprintf('  Summary for %s:\n', subjects{s});
        for b = 1:n_bands
            bn = band_names{b};
            fprintf('    %s: W=%.4f, P=%.4f, S=%.4f\n', ...
                bn, PAC_all.([bn '_words']), ...
                PAC_all.([bn '_subject_phrases']), ...
                PAC_all.([bn '_sentences']));
        end
        
        clear data word_pac_values phrase_pac_values sentence_pac_values;
        
    catch ME
        fprintf('ERROR in PAC computation for %s:\n', subjects{s});
        fprintf('  %s\n', ME.message);
        fprintf('  %s\n', ME.stack(1).name);
        continue;
    end
end

fprintf('\n=== ALL PAC COMPUTATION COMPLETE! ===\n');

%% ================= HELPER FUNCTION =================
function pac_matrix = compute_pac_all_bands(epoch_data, b_phase, a_phase, b_amp, a_amp, band_names, n_bins)
    [n_chans, ~] = size(epoch_data);
    n_bands = length(band_names);
    pac_matrix = zeros(n_bands, n_chans);
    
    phase_bins = linspace(-pi, pi, n_bins+1);
    uniform_dist = ones(1, n_bins) / n_bins;
    log_nbins = log(n_bins);
    
    for ch = 1:n_chans
        signal = epoch_data(ch, :);
        if length(signal) < 10, continue; end
        
        try
            filtered_phase = filtfilt(b_phase, a_phase, signal);
            phase_analytic = hilbert(filtered_phase);
            phase_angles = angle(phase_analytic);
        catch
            continue;
        end
        
        bin_idx = discretize(phase_angles, phase_bins);
        bin_idx(bin_idx > n_bins) = n_bins;
        
        for b = 1:n_bands
            try
                filtered_amp = filtfilt(b_amp.(band_names{b}), a_amp.(band_names{b}), signal);
                amp_envelope = abs(hilbert(filtered_amp));
            catch
                continue;
            end
            
            mean_amp = zeros(n_bins, 1);
            for bin = 1:n_bins
                idx = bin_idx == bin;
                if any(idx)
                    mean_amp(bin) = mean(amp_envelope(idx));
                end
            end
            
            if sum(mean_amp) > 0
                mean_amp = mean_amp / sum(mean_amp);
            else
                mean_amp = ones(size(mean_amp)) / n_bins;
            end
            
            mean_amp(mean_amp == 0) = eps;
            kl_div = sum(mean_amp .* log(mean_amp ./ uniform_dist'));
            pac_matrix(b, ch) = kl_div / log_nbins;
        end
    end
end