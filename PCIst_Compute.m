%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    SCRIPT 4: PCIst_Compute.m                             %%
%%    - Quantile-based discretization (3 states)            %%
%%    - 5 PCA components                                    %%
%%    - 320 ms word epochs                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0) Suppress PCA warnings
warning('off', 'stats:pca:ColRankDefX');

%% 1) Setup
clear all; 
close all;
clc;

directory = 'C:\Users\Austin\Data\';
PREPROC_Folder = 'EEG_Data\PREPROC\';
PCIst_Folder = 'PCIst\';

if ~exist([directory, PCIst_Folder], 'dir')
    mkdir([directory, PCIst_Folder]);
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

% PCIst parameters 
n_components = 5;
k_states = 3;       % quantile-based tertiles

%% 3) Subject loop
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
        
        samples_per_word = round(word_duration * Fs);
        fprintf('  Samples per word epoch: %d\n', samples_per_word);
        
        %% Compute PCIst
        fprintf('  Computing PCIst for all word epochs...\n');
        
        n_sentences_total = n_trials * num_sentences_per_trial;
        PCIst_per_sentence = zeros(n_sentences_total, words_per_sentence);
        
        sentence_counter = 0;
        
        tic;
        for tr = 1:n_trials
            trial_data = data.trial{tr};
            n_samples = size(trial_data, 2);
            
            for sent = 1:num_sentences_per_trial
                sentence_counter = sentence_counter + 1;
                
                for w = 1:words_per_sentence
                    word_idx = (sent-1) * words_per_sentence + w;
                    word_onset = trial_start_offset + (word_idx - 1) * word_duration;
                    
                    onset_sample = round(word_onset * Fs) + 1;
                    end_sample = onset_sample + samples_per_word - 1;
                    
                    if end_sample > n_samples, end_sample = n_samples; end
                    if onset_sample < 1, onset_sample = 1; end
                    
                    epoch_data = trial_data(chan_idx, onset_sample:end_sample);
                    epoch_data = epoch_data - mean(epoch_data, 2);
                    
                    PCIst_per_sentence(sentence_counter, w) = ...
                        compute_pcist(epoch_data, n_components, k_states);
                end
            end
        end
        pcist_time = toc;
        
        fprintf('    PCIst done in %.1f seconds (%d sentences)\n', ...
            pcist_time, sentence_counter);
        
        PCIst_per_sentence = PCIst_per_sentence(1:sentence_counter, :);
        
        %% Average across sentences
        PCIst_W1 = mean(PCIst_per_sentence(:, 1));
        PCIst_W2 = mean(PCIst_per_sentence(:, 2));
        PCIst_W3 = mean(PCIst_per_sentence(:, 3));
        PCIst_W4 = mean(PCIst_per_sentence(:, 4));
        
        PCIst_slope = compute_slope([PCIst_W1, PCIst_W2, PCIst_W3, PCIst_W4]);
        PCIst_W4_minus_W1 = PCIst_W4 - PCIst_W1;
        
        %% Display
        fprintf('\n  PCIst Results for %s:\n', subjects{s});
        fprintf('    Word 1: %.4f\n', PCIst_W1);
        fprintf('    Word 2: %.4f\n', PCIst_W2);
        fprintf('    Word 3: %.4f\n', PCIst_W3);
        fprintf('    Word 4: %.4f\n', PCIst_W4);
        fprintf('    Slope: %.6f\n', PCIst_slope);
        fprintf('    W4 - W1: %.4f\n', PCIst_W4_minus_W1);
        
        %% Save
        save([directory, PCIst_Folder, 'PCIst_', char(subjects(s)), '.mat'], ...
            'PCIst_W1', 'PCIst_W2', 'PCIst_W3', 'PCIst_W4', ...
            'PCIst_slope', 'PCIst_W4_minus_W1', ...
            'PCIst_per_sentence', 'Fs', ...
            'n_components', 'k_states', '-v7.3');
        
        clear data PCIst_per_sentence;
        
    catch ME
        fprintf('ERROR in PCIst computation for %s:\n', subjects{s});
        fprintf('  %s\n', ME.message);
        fprintf('  %s\n', ME.stack(1).name);
        continue;
    end
end

fprintf('\n=== ALL PCIst COMPUTATION COMPLETE! ===\n');

warning('on', 'stats:pca:ColRankDefX');

%% ================= HELPER FUNCTIONS =================

function pcist_val = compute_pcist(epoch_data, n_components, k_states)
    % Compute PCIst using quantile-based discretization 
    
    [n_chans, n_timepoints] = size(epoch_data);
    
    if n_timepoints < 10 || n_chans < 2
        pcist_val = 0;
        return;
    end
    
    try
        X = epoch_data';
        [~, score, ~, ~, ~] = pca(X);
        
        n_comp = min([n_components, size(score, 2)]);
        components = score(:, 1:n_comp);
        
        total_transitions = 0;
        
        for c = 1:n_comp
            sig = components(:, c);
            
            if std(sig) <= 0
                continue;
            end
            
            % Quantile-based binning into k_states
            edges = quantile(sig, linspace(0, 1, k_states + 1));
            edges(1) = -inf;
            edges(end) = inf;
            
            states = discretize(sig, edges);
            
            valid = ~isnan(states);
            states = states(valid);
            
            if length(states) < 2
                continue;
            end
            
            transitions = sum(states(2:end) ~= states(1:end-1));
            total_transitions = total_transitions + transitions;
        end
        
        if n_comp > 0 && n_timepoints > 1
            pcist_val = total_transitions / (n_comp * (n_timepoints - 1));
        else
            pcist_val = 0;
        end
        
    catch
        pcist_val = 0;
    end
end

function slope = compute_slope(values)
    x = 1:length(values);
    p = polyfit(x, values, 1);
    slope = p(1);
end