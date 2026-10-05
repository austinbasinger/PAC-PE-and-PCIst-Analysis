%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    SCRIPT 3: Permutation Entropy (PE) Computation        %%
%%    - Reads data.trial{tr} (channels x time)              %%
%%    - Resamples to 256 Hz if needed                       %%
%%    - PE per channel, averaged across channels            %%
%%    - 4 values per patient (Word 1, 2, 3, 4)              %%
%%    - Gradient metrics (slope, W4-W1)                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1) Define data directories, subjects and add paths
clear all; 
close all;
clc;

directory = 'C:\Users\Austin\Data\';
PREPROC_Folder = 'EEG_Data\PREPROC\';
PE_Folder = 'PE\';

if ~exist([directory, PE_Folder], 'dir')
    mkdir([directory, PE_Folder]);
end

% Sokoliuk's preprocessed files are P01-P17
subjects = [{'P01'},{'P02'},{'P03'},{'P04'},{'P05'},{'P06'},{'P07'},{'P08'},...
    {'P09'},{'P10'},{'P11'},{'P12'},{'P13'},{'P14'},{'P15'},{'P16'},{'P17'}];

% Map to original raw labels
raw_labels = [{'P04'},{'P06'},{'P08'},{'P10'},{'P12'},{'P13'},{'P14'},{'P15'},...
    {'P17'},{'P20'},{'P21'},{'P22'},{'P24'},{'P25'},{'P26'},{'P27'},{'P28'}];

%% 2) Define stimulus parameters
word_duration = 0.32;
words_per_sentence = 4;
num_sentences_per_trial = 12;
trial_start_offset = 1.0;
chans2exclude = {'A1', 'A2', 'COMNT', 'SCALE'};
TARGET_FS = 256;

% Permutation Entropy parameters
m = 3;      % Embedding dimension (order of patterns)
tau = 1;    % Time delay between samples

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
        
        %% 3a) Resample to TARGET_FS if needed
        if Fs_original ~= TARGET_FS
            fprintf('  Resampling from %d Hz to %d Hz...\n', Fs_original, TARGET_FS);
            
            for tr = 1:n_trials
                trial_data = data.trial{tr};
                resampled = resample(trial_data', TARGET_FS, Fs_original)';
                data.trial{tr} = resampled;
            end
            
            data.fsample = TARGET_FS;
            Fs = TARGET_FS;
        else
            Fs = Fs_original;
        end
        
        fprintf('  Final sampling rate: %d Hz\n', Fs);
        
        chan_idx = find(~ismember(data.label, chans2exclude));
        n_chans = length(chan_idx);
        fprintf('  Using %d EEG channels\n', n_chans);
        
        samples_per_word = round(word_duration * Fs);
        fprintf('  Samples per word epoch: %d\n', samples_per_word);
        
        %% 3b) Compute PE for each word epoch
        fprintf('  Computing Permutation Entropy for all word epochs...\n');
        
        n_sentences_total = n_trials * num_sentences_per_trial;
        PE_per_sentence = zeros(n_sentences_total, words_per_sentence);
        
        sentence_counter = 0;
        
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
                    
                    if end_sample > n_samples
                        end_sample = n_samples;
                    end
                    if onset_sample < 1
                        onset_sample = 1;
                    end
                    
                    epoch_data = trial_data(chan_idx, onset_sample:end_sample);
                    
                    % Compute PE per channel, average across channels
                    pe_vals = zeros(n_chans, 1);
                    for ch = 1:n_chans
                        pe_vals(ch) = compute_permutation_entropy(epoch_data(ch, :), m, tau);
                    end
                    
                    PE_per_sentence(sentence_counter, w) = mean(pe_vals);
                end
            end
        end
        
        PE_per_sentence = PE_per_sentence(1:sentence_counter, :);
        fprintf('  Computed PE for %d sentences\n', sentence_counter);
        
        %% 3c) Average across sentences for each word position
        PE_W1 = mean(PE_per_sentence(:, 1));
        PE_W2 = mean(PE_per_sentence(:, 2));
        PE_W3 = mean(PE_per_sentence(:, 3));
        PE_W4 = mean(PE_per_sentence(:, 4));
        
        PE_slope = compute_slope([PE_W1, PE_W2, PE_W3, PE_W4]);
        PE_W4_minus_W1 = PE_W4 - PE_W1;
        
        %% 3d) Display results
        fprintf('\n  PE Results for %s:\n', subjects{s});
        fprintf('    Word 1: %.4f\n', PE_W1);
        fprintf('    Word 2: %.4f\n', PE_W2);
        fprintf('    Word 3: %.4f\n', PE_W3);
        fprintf('    Word 4: %.4f\n', PE_W4);
        fprintf('    Slope: %.6f\n', PE_slope);
        fprintf('    W4 - W1: %.4f\n', PE_W4_minus_W1);
        
        %% 3e) Save results
        save([directory, PE_Folder, 'PE_', char(subjects(s)), '.mat'], ...
            'PE_W1', 'PE_W2', 'PE_W3', 'PE_W4', ...
            'PE_slope', 'PE_W4_minus_W1', ...
            'PE_per_sentence', 'Fs', 'm', 'tau', '-v7.3');
        
        clear data;
        
    catch ME
        fprintf('ERROR in PE computation for %s:\n', subjects{s});
        fprintf('  %s\n', ME.message);
        continue;
    end
end

fprintf('\n=== ALL PERMUTATION ENTROPY COMPUTATION COMPLETE! ===\n');

%% 4) Helper function: Permutation Entropy
function pe = compute_permutation_entropy(signal, m, tau)
    % compute_permutation_entropy: Calculate Permutation Entropy
    % Inputs:
    %   signal: 1D signal vector
    %   m: embedding dimension (pattern length)
    %   tau: time delay between samples
    % Output:
    %   pe: Permutation Entropy (normalized between 0 and 1)
    
    N = length(signal);
    
    % Check if signal is long enough
    if N < (m - 1) * tau + 1
        pe = 0;
        return;
    end
    
    % Number of possible patterns = m!
    n_patterns = factorial(m);
    
    % Number of embedding vectors
    n_vectors = N - (m - 1) * tau;
    
    % Build embedding matrix (each row is an embedding vector)
    % Use a moving window approach
    embedded = zeros(n_vectors, m);
    for i = 1:n_vectors
        for j = 1:m
            embedded(i, j) = signal(i + (j - 1) * tau);
        end
    end
    
    % For each embedding vector, determine the ordinal pattern
    % The pattern is defined by the relative order of the values
    pattern_counts = zeros(n_patterns, 1);
    
    for i = 1:n_vectors
        vec = embedded(i, :);
        
        % Get the ranking (ordinal pattern) of this vector
        [~, rank_order] = sort(vec);
        
        % Convert rank order to a pattern index (1 to m!)
        % We use the Lehmer code (factoradic) to index patterns
        pattern_idx = ranking_to_index(rank_order, m);
        pattern_counts(pattern_idx) = pattern_counts(pattern_idx) + 1;
    end
    
    % Convert counts to probabilities
    probabilities = pattern_counts / sum(pattern_counts);
    
    % Compute Shannon entropy
    % Avoid log(0) by only summing over non-zero probabilities
    pe = 0;
    for i = 1:length(probabilities)
        if probabilities(i) > 0
            pe = pe - probabilities(i) * log2(probabilities(i));
        end
    end
    
    % Normalize by log2(m!) to get value between 0 and 1
    pe = pe / log2(n_patterns);
end

%% 5) Helper function: Convert ranking to pattern index
function idx = ranking_to_index(rank_order, m)
    % Convert a ranking vector to a unique index (1 to m!)
    % Uses Lehmer code (factoradic)
    
    available = 1:m;
    idx = 1;
    
    for i = 1:m
        % Position of rank_order(i) in the available list
        pos = find(available == rank_order(i));
        idx = idx + (pos - 1) * factorial(m - i);
        available(pos) = [];  % Remove this element
    end
end

%% 6) Helper function: Compute slope
function slope = compute_slope(values)
    x = 1:length(values);
    p = polyfit(x, values, 1);
    slope = p(1);
end