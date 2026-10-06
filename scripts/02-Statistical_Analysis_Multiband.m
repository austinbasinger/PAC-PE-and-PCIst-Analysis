%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    SCRIPT 2: Statistical_Analysis_Multiband.m            %%
%%    Multi-Band PAC Statistics on Sokoliuk's data          %%
%%    - Correlations: 3M and 6M for ALL bands × levels      %%
%%    - Bootstrap specificity tests                         %%
%%    - Backward elimination linear regression              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0) Setup
clear all; 
close all;
clc;

directory = 'C:\Users\Austin\Data\';
PAC_Folder = 'PAC_MultiBand_All\';
FIGURES_Folder = 'Figures_PAC_MultiBand\';
RESULTS_Folder = 'Results_PAC_MultiBand\';

if ~exist([directory, FIGURES_Folder], 'dir')
    mkdir([directory, FIGURES_Folder]);
end
if ~exist([directory, RESULTS_Folder], 'dir')
    mkdir([directory, RESULTS_Folder]);
end

% Sokoliuk's preprocessed subject labels
subjects = [{'P01'},{'P02'},{'P03'},{'P04'},{'P05'},{'P06'},{'P07'},{'P08'},...
    {'P09'},{'P10'},{'P11'},{'P12'},{'P13'},{'P14'},{'P15'},{'P16'},{'P17'}];

raw_labels = [{'P04'},{'P06'},{'P08'},{'P10'},{'P12'},{'P13'},{'P14'},{'P15'},...
    {'P17'},{'P20'},{'P21'},{'P22'},{'P24'},{'P25'},{'P26'},{'P27'},{'P28'}];

% GOSE outcomes mapped to P01-P17
GOSE_3M = [1 1 2 3 3 3 2 1 3 6 3 3 6 3 3 4 4];
GOSE_6M = [1 1 2 3 3 3 3 1 2 7 3 3 6 NaN 3 6 5];

% Clinical variables mapped to P01-P17
clinical_data.age = [72 86 26 40 59 44 82 64 70 70 27 77 54 59 59 61 32];
clinical_data.gcs = [3 4 4 3 1 4 1 3 4 5 4 4 4 4 3 3 5];
clinical_data.days = [5 5 17 12 13 10 3 9 5 10 19 12 10 9 14 15 17];
clinical_data.ct = [2 2 5 5 5 5 5 5 5 6 2 2 2 3 5 2 2];

n_subjects = length(subjects);
subs_3M = 1:n_subjects;
subs_6M = find(~isnan(GOSE_6M));

fprintf('=== PAC MULTIBAND STATISTICAL ANALYSIS ===\n');
fprintf('Subjects: %d (6M: %d)\n\n', n_subjects, length(subs_6M));

%% 1) Define bands and levels
bands = {'delta', 'theta', 'alpha', 'beta', 'low_gamma', 'mid_gamma', 'high_gamma'};
band_labels = {'Delta (1-4 Hz)', 'Theta (4-8 Hz)', 'Alpha (8-13 Hz)', ...
    'Beta (13-30 Hz)', 'Low Gamma (30-45 Hz)', 'Mid Gamma (45-60 Hz)', ...
    'High Gamma (65-80 Hz)'};
levels = {'words', 'subject_phrases', 'sentences'};
level_labels = {'Words', 'Subject Phrases', 'Sentences'};

%% 2) Load PAC data
fprintf('Loading PAC data...\n');

PAC_data = zeros(n_subjects, length(bands), length(levels));

for s = 1:n_subjects
    f = [directory, PAC_Folder, 'PAC_MultiBand_All_', subjects{s}, '.mat'];
    if exist(f, 'file')
        load(f);
        for b = 1:length(bands)
            for l = 1:length(levels)
                field_name = [bands{b} '_' levels{l}];
                if isfield(PAC_all, field_name)
                    PAC_data(s, b, l) = PAC_all.(field_name);
                else
                    PAC_data(s, b, l) = NaN;
                end
            end
        end
    else
        warning('PAC file not found for %s', subjects{s});
        PAC_data(s, :, :) = NaN;
    end
end

fprintf('Loaded PAC data for %d subjects\n\n', sum(~isnan(PAC_data(:,1,1))));

%% 3) Correlations - All bands × All levels
fprintf('============================================\n');
fprintf('=== CORRELATIONS ===\n');
fprintf('============================================\n\n');

fprintf('%-20s %-20s %-15s %-15s %-15s %-15s\n', ...
    'Band', 'Level', '3M rho', '3M p', '6M rho', '6M p');
fprintf('%-20s %-20s %-15s %-15s %-15s %-15s\n', ...
    '--------------------', '--------------------', '---------------', '---------------', '---------------', '---------------');

corr_results = struct();

for b = 1:length(bands)
    for l = 1:length(levels)
        pac_vals = squeeze(PAC_data(:, b, l));
        
        [rho_3M, p_3M] = corr(GOSE_3M(subs_3M)', pac_vals(subs_3M), ...
            'Type', 'Spearman', 'Rows', 'complete');
        [rho_6M, p_6M] = corr(GOSE_6M(subs_6M)', pac_vals(subs_6M), ...
            'Type', 'Spearman', 'Rows', 'complete');
        
        corr_results.([bands{b} '_' levels{l} '_3M']) = struct('rho', rho_3M, 'p', p_3M);
        corr_results.([bands{b} '_' levels{l} '_6M']) = struct('rho', rho_6M, 'p', p_6M);
        
        sig_3M = ''; sig_6M = '';
        if p_3M < 0.05, sig_3M = '*'; end
        if p_6M < 0.05, sig_6M = '*'; end
        
        fprintf('%-20s %-20s %-15.3f%s %-15.4f %-15.3f%s %-15.4f\n', ...
            band_labels{b}, level_labels{l}, rho_3M, sig_3M, p_3M, rho_6M, sig_6M, p_6M);
    end
end

fprintf('\n* p < 0.05\n');

%% 4) Correlation Figures for each band × level
fprintf('\n=== Creating Correlation Figures ===\n');

for b = 1:length(bands)
    for l = 1:length(levels)
        pac_vals = squeeze(PAC_data(:, b, l));
        
        [rho_3M, p_3M] = corr(GOSE_3M(subs_3M)', pac_vals(subs_3M), ...
            'Type', 'Spearman', 'Rows', 'complete');
        [rho_6M, p_6M] = corr(GOSE_6M(subs_6M)', pac_vals(subs_6M), ...
            'Type', 'Spearman', 'Rows', 'complete');
        
        % 3M
        figure('Position', [100 100 600 500]);
        hold on;
        scatter(GOSE_3M(subs_3M), pac_vals(subs_3M), 120, 'b', 'filled', 'MarkerEdgeColor', 'k');
        p_fit = polyfit(GOSE_3M(subs_3M), pac_vals(subs_3M), 1);
        x_vals = linspace(min(GOSE_3M(subs_3M)), max(GOSE_3M(subs_3M)), 100);
        plot(x_vals, polyval(p_fit, x_vals), 'r-', 'LineWidth', 2);
        xlabel('GOSE 3 Months', 'FontSize', 14);
        ylabel(sprintf('%s PAC - %s', band_labels{b}, level_labels{l}), 'FontSize', 12);
        title(sprintf('%s %s PAC vs GOSE 3M\nrho = %.3f, p = %.4f', ...
            band_labels{b}, level_labels{l}, rho_3M, p_3M), 'FontSize', 12);
        xlim([0 8]); grid on; set(gca, 'XTick', 1:8);
        saveas(gcf, [directory, FIGURES_Folder, sprintf('PAC_%s_%s_3M.png', bands{b}, levels{l})], 'png');
        close(gcf);
        
        % 6M
        figure('Position', [100 100 600 500]);
        hold on;
        scatter(GOSE_6M(subs_6M), pac_vals(subs_6M), 120, 'b', 'filled', 'MarkerEdgeColor', 'k');
        p_fit = polyfit(GOSE_6M(subs_6M), pac_vals(subs_6M), 1);
        x_vals = linspace(min(GOSE_6M(subs_6M)), max(GOSE_6M(subs_6M)), 100);
        plot(x_vals, polyval(p_fit, x_vals), 'r-', 'LineWidth', 2);
        xlabel('GOSE 6 Months', 'FontSize', 14);
        ylabel(sprintf('%s PAC - %s', band_labels{b}, level_labels{l}), 'FontSize', 12);
        title(sprintf('%s %s PAC vs GOSE 6M\nrho = %.3f, p = %.4f', ...
            band_labels{b}, level_labels{l}, rho_6M, p_6M), 'FontSize', 12);
        xlim([0 8]); grid on; set(gca, 'XTick', 1:8);
        saveas(gcf, [directory, FIGURES_Folder, sprintf('PAC_%s_%s_6M.png', bands{b}, levels{l})], 'png');
        close(gcf);
    end
end

fprintf('Figures saved.\n');

%% 5) Bootstrap for Word-level PAC of each band
fprintf('\n=== Bootstrap Tests (Word PAC) ===\n\n');

n_reps = 1000;
bootstrap_results = struct();

for b = 1:length(bands)
    pac_vals = squeeze(PAC_data(:, b, 1));  % Word level
    
    [rho_3M, ~] = corr(GOSE_3M(subs_3M)', pac_vals(subs_3M), 'Type', 'Spearman');
    [rho_6M, ~] = corr(GOSE_6M(subs_6M)', pac_vals(subs_6M), 'Type', 'Spearman');
    
    surr_3M = zeros(n_reps, 1);
    for rep = 1:n_reps
        shuffled = pac_vals(randperm(n_subjects));
        [surr_3M(rep), ~] = corr(GOSE_3M(subs_3M)', shuffled(subs_3M), 'Type', 'Spearman');
    end
    p_boot_3M = sum(abs(surr_3M) >= abs(rho_3M)) / n_reps;
    
    surr_6M = zeros(n_reps, 1);
    for rep = 1:n_reps
        shuffled = pac_vals(randperm(n_subjects));
        [surr_6M(rep), ~] = corr(GOSE_6M(subs_6M)', shuffled(subs_6M), 'Type', 'Spearman');
    end
    p_boot_6M = sum(abs(surr_6M) >= abs(rho_6M)) / n_reps;
    
    bootstrap_results.(bands{b}) = struct('p_boot_3M', p_boot_3M, 'p_boot_6M', p_boot_6M);
    
    fprintf('%-20s : 3M p_boot = %.4f, 6M p_boot = %.4f\n', ...
        band_labels{b}, p_boot_3M, p_boot_6M);
end

%% 6) Backward Elimination Regression - Word-level PAC for each band
fprintf('\n============================================\n');
fprintf('=== BACKWARD ELIMINATION REGRESSION ===\n');
fprintf('============================================\n\n');

reg_results = struct();

for t = 1:2
    if t == 1
        timepoint = '3M';
        subs = subs_3M;
        Y = GOSE_3M(subs)';
    else
        timepoint = '6M';
        subs = subs_6M;
        Y = GOSE_6M(subs)';
    end
    n_reg = length(subs);
    
    fprintf('\n--- %s OUTCOME (n = %d) ---\n', timepoint, n_reg);
    
    X_age = clinical_data.age(subs)';
    X_gcs = clinical_data.gcs(subs)';
    X_days = clinical_data.days(subs)';
    X_ct = clinical_data.ct(subs)';
    
    Y_norm = norminv((tiedrank(Y) - 0.5) / length(Y));
    
    for b = 1:length(bands)
        pac_vals = squeeze(PAC_data(:, b, 1));  % Word level
        
        X_full = [ones(n_reg,1), X_age, X_gcs, X_days, X_ct, pac_vals(subs)];
        pred_names = {'(Intercept)', 'Age', 'GCS', 'Days', 'CT', ...
            sprintf('PAC_%s_words', bands{b})};
        
        X_curr = X_full;
        pred_curr = pred_names;
        removed = {};
        final_R2 = 0; final_adjR2 = 0; final_F = NaN; final_p = NaN;
        final_preds = {};
        
        while true
            [bb, ~, ~, ~, stats] = regress(Y_norm, X_curr);
            R2 = stats(1); F = stats(2); p_model = stats(3);
            adjR2 = 1 - (1-R2) * (n_reg-1) / (n_reg - size(X_curr,2));
            
            resid = Y_norm - X_curr*bb;
            s2 = sum(resid.^2) / (n_reg - size(X_curr,2));
            cov_b = s2 * inv(X_curr'*X_curr);
            se_b = sqrt(diag(cov_b));
            t_vals = bb ./ se_b;
            p_vals = 2 * (1 - tcdf(abs(t_vals), n_reg - size(X_curr,2)));
            
            if length(pred_curr) <= 1
                final_R2 = R2; final_adjR2 = adjR2;
                final_F = F; final_p = p_model;
                final_preds = {};
                break;
            end
            
            p_nonint = p_vals(2:end);
            [max_p, max_idx] = max(p_nonint);
            
            if max_p > 0.1
                remove_idx = max_idx + 1;
                removed{end+1} = pred_curr{remove_idx};
                keep = setdiff(1:length(pred_curr), remove_idx);
                X_curr = X_curr(:, keep);
                pred_curr = pred_curr(keep);
            else
                final_R2 = R2; final_adjR2 = adjR2;
                final_F = F; final_p = p_model;
                final_preds = pred_curr(2:end);
                break;
            end
        end
        
        reg_results.([bands{b} '_' timepoint]) = struct(...
            'R2', final_R2, 'adjR2', final_adjR2, ...
            'F', final_F, 'p', final_p, ...
            'predictors', {final_preds}, 'removed', {removed});
        
        pred_str = strjoin(final_preds, ', ');
        if isempty(pred_str), pred_str = '(none)'; end
        
        fprintf('  %s: R2 = %.3f, adjR2 = %.3f, p = %.4f | Retained: %s\n', ...
            band_labels{b}, final_R2, final_adjR2, final_p, pred_str);
    end
end

%% 7) Save
save([directory, RESULTS_Folder, 'PAC_MultiBand_Statistics.mat'], ...
    'corr_results', 'bootstrap_results', 'reg_results', 'PAC_data', ...
    'GOSE_3M', 'GOSE_6M', 'subjects', 'raw_labels');

fprintf('\n=== PAC ANALYSIS COMPLETE ===\n');
fprintf('Results saved in: %s%s\n', directory, RESULTS_Folder);
fprintf('Figures saved in: %s%s\n', directory, FIGURES_Folder);
