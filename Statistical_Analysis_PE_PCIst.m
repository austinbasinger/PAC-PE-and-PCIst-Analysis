%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    SCRIPT 5: Statistical Analysis - PE and PCIst         %%
%%    - Correlations: slope and W4-W1 with GOSE 3M and 6M   %%
%%    - Bootstrap specificity for each metric               %%
%%    - Backward elimination regression                     %%
%%    - Comparison: PE vs PCIst                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0) Setup
clear all; 
close all;
clc;

warning('off', 'stats:pca:ColRankDefX');

directory = 'C:\Users\Austin\Data\';
PE_Folder = 'PE\';
PCIst_Folder = 'PCIst\';
FIGURES_Folder = 'Figures_PE_PCIst\';
RESULTS_Folder = 'Results_PE_PCIst\';

if ~exist([directory, FIGURES_Folder], 'dir')
    mkdir([directory, FIGURES_Folder]);
end
if ~exist([directory, RESULTS_Folder], 'dir')
    mkdir([directory, RESULTS_Folder]);
end

% Preprocessed subject labels (P01-P17)
subjects = [{'P01'},{'P02'},{'P03'},{'P04'},{'P05'},{'P06'},{'P07'},{'P08'},...
    {'P09'},{'P10'},{'P11'},{'P12'},{'P13'},{'P14'},{'P15'},{'P16'},{'P17'}];

% Raw labels for reference
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

fprintf('=== STATISTICAL ANALYSIS: PE vs PCIst ===\n');
fprintf('Total subjects: %d (6-month: %d)\n\n', n_subjects, length(subs_6M));

%% 1) Load PE and PCIst data
fprintf('Loading PE and PCIst data...\n');

PE = struct();
PCIst = struct();

for s = 1:n_subjects
    % PE
    pe_file = [directory, PE_Folder, 'PE_', subjects{s}, '.mat'];
    if exist(pe_file, 'file')
        load(pe_file);
        PE.W1(s) = PE_W1;
        PE.W2(s) = PE_W2;
        PE.W3(s) = PE_W3;
        PE.W4(s) = PE_W4;
        PE.slope(s) = PE_slope;
        PE.W4_minus_W1(s) = PE_W4_minus_W1;
    else
        PE.W1(s) = NaN; PE.W2(s) = NaN; PE.W3(s) = NaN; PE.W4(s) = NaN;
        PE.slope(s) = NaN; PE.W4_minus_W1(s) = NaN;
    end
    
    % PCIst
    pcist_file = [directory, PCIst_Folder, 'PCIst_', subjects{s}, '.mat'];
    if exist(pcist_file, 'file')
        load(pcist_file);
        PCIst.W1(s) = PCIst_W1;
        PCIst.W2(s) = PCIst_W2;
        PCIst.W3(s) = PCIst_W3;
        PCIst.W4(s) = PCIst_W4;
        PCIst.slope(s) = PCIst_slope;
        PCIst.W4_minus_W1(s) = PCIst_W4_minus_W1;
    else
        PCIst.W1(s) = NaN; PCIst.W2(s) = NaN; PCIst.W3(s) = NaN; PCIst.W4(s) = NaN;
        PCIst.slope(s) = NaN; PCIst.W4_minus_W1(s) = NaN;
    end
end

fprintf('PE loaded for %d subjects\n', sum(~isnan(PE.W1)));
fprintf('PCIst loaded for %d subjects\n\n', sum(~isnan(PCIst.W1)));

%% 2) Correlations: slope and W4-W1 with GOSE 3M and 6M
fprintf('============================================\n');
fprintf('=== CORRELATIONS ===\n');
fprintf('============================================\n\n');

corr_results = struct();

measures = {'PE', 'PCIst'};
metrics = {'slope', 'W4_minus_W1'};

fprintf('%-8s %-12s %-15s %-15s %-15s %-15s\n', ...
    'Measure', 'Metric', '3M rho', '3M p', '6M rho', '6M p');
fprintf('%-8s %-12s %-15s %-15s %-15s %-15s\n', ...
    '--------', '----------', '---------------', '---------------', '---------------', '---------------');

for m = 1:length(measures)
    measure = measures{m};
    if strcmp(measure, 'PE')
        data = PE;
    else
        data = PCIst;
    end
    
    for k = 1:length(metrics)
        metric = metrics{k};
        vals = data.(metric);
        
        [rho_3M, p_3M] = corr(GOSE_3M(subs_3M)', vals(subs_3M)', ...
            'Type', 'Spearman', 'Rows', 'complete');
        
        [rho_6M, p_6M] = corr(GOSE_6M(subs_6M)', vals(subs_6M)', ...
            'Type', 'Spearman', 'Rows', 'complete');
        
        corr_results.([measure '_' metric '_3M']) = struct('rho', rho_3M, 'p', p_3M);
        corr_results.([measure '_' metric '_6M']) = struct('rho', rho_6M, 'p', p_6M);
        
        sig_3M = ''; sig_6M = '';
        if p_3M < 0.05, sig_3M = '*'; end
        if p_6M < 0.05, sig_6M = '*'; end
        
        fprintf('%-8s %-12s %-15.3f%s %-15.4f %-15.3f%s %-15.4f\n', ...
            measure, strrep(metric,'_','-'), rho_3M, sig_3M, p_3M, rho_6M, sig_6M, p_6M);
    end
end

fprintf('\n* p < 0.05\n');

%% 3) Correlation Figures
fprintf('\n=== Creating Correlation Figures ===\n');

for m = 1:length(measures)
    measure = measures{m};
    if strcmp(measure, 'PE')
        data = PE;
    else
        data = PCIst;
    end
    
    for k = 1:length(metrics)
        metric = metrics{k};
        vals = data.(metric);
        metric_label = strrep(metric, '_', '-');
        
        % ---- 3M Figure ----
        [rho_3M, p_3M] = corr(GOSE_3M(subs_3M)', vals(subs_3M)', ...
            'Type', 'Spearman', 'Rows', 'complete');
        
        figure('Position', [100 100 600 500]);
        hold on;
        scatter(GOSE_3M(subs_3M), vals(subs_3M), 120, 'b', 'filled', ...
            'MarkerEdgeColor', 'k');
        p_fit = polyfit(GOSE_3M(subs_3M), vals(subs_3M), 1);
        x_vals = linspace(min(GOSE_3M(subs_3M)), max(GOSE_3M(subs_3M)), 100);
        plot(x_vals, polyval(p_fit, x_vals), 'r-', 'LineWidth', 2);
        xlabel('GOSE 3 Months', 'FontSize', 14);
        ylabel(sprintf('%s %s', measure, metric_label), 'FontSize', 14);
        title(sprintf('%s %s vs GOSE 3M\nrho = %.3f, p = %.4f', ...
            measure, metric_label, rho_3M, p_3M), 'FontSize', 12);
        xlim([0 8]);
        grid on;
        set(gca, 'XTick', 1:8);
        saveas(gcf, [directory, FIGURES_Folder, ...
            sprintf('%s_%s_GOSE_3M.png', measure, metric)], 'png');
        close(gcf);
        
        % ---- 6M Figure ----
        [rho_6M, p_6M] = corr(GOSE_6M(subs_6M)', vals(subs_6M)', ...
            'Type', 'Spearman', 'Rows', 'complete');
        
        figure('Position', [100 100 600 500]);
        hold on;
        scatter(GOSE_6M(subs_6M), vals(subs_6M), 120, 'b', 'filled', ...
            'MarkerEdgeColor', 'k');
        p_fit = polyfit(GOSE_6M(subs_6M), vals(subs_6M), 1);
        x_vals = linspace(min(GOSE_6M(subs_6M)), max(GOSE_6M(subs_6M)), 100);
        plot(x_vals, polyval(p_fit, x_vals), 'r-', 'LineWidth', 2);
        xlabel('GOSE 6 Months', 'FontSize', 14);
        ylabel(sprintf('%s %s', measure, metric_label), 'FontSize', 14);
        title(sprintf('%s %s vs GOSE 6M\nrho = %.3f, p = %.4f', ...
            measure, metric_label, rho_6M, p_6M), 'FontSize', 12);
        xlim([0 8]);
        grid on;
        set(gca, 'XTick', 1:8);
        
        for i = 1:length(subs_6M)
            idx = subs_6M(i);
            text(GOSE_6M(idx)+0.1, vals(idx), subjects{idx}, 'FontSize', 8);
        end
        
        saveas(gcf, [directory, FIGURES_Folder, ...
            sprintf('%s_%s_GOSE_6M.png', measure, metric)], 'png');
        close(gcf);
    end
end

fprintf('Figures saved.\n');

%% 4) Bootstrap Specificity for Each Metric
fprintf('\n=== Bootstrap Specificity Tests ===\n');
fprintf('(1000 iterations per metric)\n\n');

n_reps = 1000;
bootstrap_results = struct();

for m = 1:length(measures)
    measure = measures{m};
    if strcmp(measure, 'PE')
        data = PE;
    else
        data = PCIst;
    end
    
    for k = 1:length(metrics)
        metric = metrics{k};
        vals = data.(metric);
        
        % Get actual correlations
        [rho_3M, ~] = corr(GOSE_3M(subs_3M)', vals(subs_3M)', 'Type', 'Spearman');
        [rho_6M, ~] = corr(GOSE_6M(subs_6M)', vals(subs_6M)', 'Type', 'Spearman');
        
        % Bootstrap 3M - shuffle ALL subjects, then subset
        surr_3M = zeros(n_reps, 1);
        for rep = 1:n_reps
            shuffled = vals(randperm(n_subjects));
            [surr_3M(rep), ~] = corr(GOSE_3M(subs_3M)', shuffled(subs_3M)', ...
                'Type', 'Spearman');
        end
        p_boot_3M = sum(abs(surr_3M) >= abs(rho_3M)) / n_reps;
        
        % Bootstrap 6M - shuffle ALL subjects, then subset
        surr_6M = zeros(n_reps, 1);
        for rep = 1:n_reps
            shuffled = vals(randperm(n_subjects));
            [surr_6M(rep), ~] = corr(GOSE_6M(subs_6M)', shuffled(subs_6M)', ...
                'Type', 'Spearman');
        end
        p_boot_6M = sum(abs(surr_6M) >= abs(rho_6M)) / n_reps;
        
        bootstrap_results.([measure '_' metric '_3M']) = p_boot_3M;
        bootstrap_results.([measure '_' metric '_6M']) = p_boot_6M;
        
        fprintf('%-8s %-12s : 3M p_boot = %.4f, 6M p_boot = %.4f\n', ...
            measure, strrep(metric,'_','-'), p_boot_3M, p_boot_6M);
        
        % Plot bootstrap for 6M
        figure('Position', [100 100 700 500]);
        histogram(surr_6M, 30, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k');
        hold on;
        line([rho_6M rho_6M], ylim, 'Color', 'r', 'LineWidth', 3);
        xlabel('Surrogate rho', 'FontSize', 14);
        ylabel('Frequency', 'FontSize', 14);
        title(sprintf('Bootstrap: %s %s vs GOSE 6M (p = %.4f)', ...
            measure, strrep(metric,'_','-'), p_boot_6M), 'FontSize', 12);
        legend({'Surrogate', 'Actual'}, 'FontSize', 12);
        grid on;
        saveas(gcf, [directory, FIGURES_Folder, ...
            sprintf('Bootstrap_%s_%s_6M.png', measure, metric)], 'png');
        close(gcf);
        
        % Plot bootstrap for 3M
        figure('Position', [100 100 700 500]);
        histogram(surr_3M, 30, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k');
        hold on;
        line([rho_3M rho_3M], ylim, 'Color', 'r', 'LineWidth', 3);
        xlabel('Surrogate rho', 'FontSize', 14);
        ylabel('Frequency', 'FontSize', 14);
        title(sprintf('Bootstrap: %s %s vs GOSE 3M (p = %.4f)', ...
            measure, strrep(metric,'_','-'), p_boot_3M), 'FontSize', 12);
        legend({'Surrogate', 'Actual'}, 'FontSize', 12);
        grid on;
        saveas(gcf, [directory, FIGURES_Folder, ...
            sprintf('Bootstrap_%s_%s_3M.png', measure, metric)], 'png');
        close(gcf);
    end
end

%% 5) Backward Elimination Regression (with diagnostics)
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
    
    for m = 1:length(measures)
        measure = measures{m};
        if strcmp(measure, 'PE')
            data = PE;
        else
            data = PCIst;
        end
        
        for k = 1:length(metrics)
            metric = metrics{k};
            vals = data.(metric);
            
            X_full = [ones(n_reg,1), X_age, X_gcs, X_days, X_ct, vals(subs)'];
            pred_names = {'(Intercept)', 'Age', 'GCS', 'Days', 'CT', ...
                sprintf('%s_%s', measure, metric)};
            
            X_curr = X_full;
            pred_curr = pred_names;
            removed = {};
            final_R2 = 0;
            final_adjR2 = 0;
            final_F = NaN;
            final_p = NaN;
            final_preds = {};
            step = 0;
            
            fprintf('\n  --- %s %s ---\n', measure, strrep(metric,'_','-'));
            
            while true
                step = step + 1;
                [b, ~, ~, ~, stats] = regress(Y_norm, X_curr);
                R2 = stats(1); F = stats(2); p_model = stats(3);
                adjR2 = 1 - (1-R2) * (n_reg-1) / (n_reg - size(X_curr,2));
                
                resid = Y_norm - X_curr*b;
                s2 = sum(resid.^2) / (n_reg - size(X_curr,2));
                cov_b = s2 * inv(X_curr'*X_curr);
                se_b = sqrt(diag(cov_b));
                t_vals = b ./ se_b;
                p_vals = 2 * (1 - tcdf(abs(t_vals), n_reg - size(X_curr,2)));
                
                % Diagnostic output
                fprintf('    Step %d: R2 = %.3f, adjR2 = %.3f, F = %.3f, p = %.4f\n', ...
                    step, R2, adjR2, F, p_model);
                fprintf('      Predictors and p-values:\n');
                for i = 1:length(pred_curr)
                    fprintf('        %-25s beta = %8.4f, p = %.4f\n', ...
                        pred_curr{i}, b(i), p_vals(i));
                end
                
                if length(pred_curr) <= 1
                    fprintf('      -> Only intercept remains. Stopping.\n');
                    final_R2 = R2; final_adjR2 = adjR2;
                    final_F = F; final_p = p_model;
                    final_preds = {};
                    break;
                end
                
                p_nonint = p_vals(2:end);
                [max_p, max_idx] = max(p_nonint);
                
                if max_p > 0.1
                    remove_idx = max_idx + 1;
                    fprintf('      -> Removing "%s" (p = %.4f > 0.1)\n', ...
                        pred_curr{remove_idx}, max_p);
                    removed{end+1} = pred_curr{remove_idx};
                    keep = setdiff(1:length(pred_curr), remove_idx);
                    X_curr = X_curr(:, keep);
                    pred_curr = pred_curr(keep);
                else
                    fprintf('      -> All remaining predictors have p < 0.1. Stopping.\n');
                    final_R2 = R2; final_adjR2 = adjR2;
                    final_F = F; final_p = p_model;
                    final_preds = pred_curr(2:end);
                    break;
                end
            end
            
            reg_results.([measure '_' metric '_' timepoint]) = struct(...
                'R2', final_R2, 'adjR2', final_adjR2, ...
                'F', final_F, 'p', final_p, ...
                'predictors', {final_preds}, 'removed', {removed});
            
            pred_str = strjoin(final_preds, ', ');
            if isempty(pred_str), pred_str = '(none)'; end
            
            fprintf('    FINAL: R2 = %.3f, adjR2 = %.3f, p = %.4f | Retained: %s\n', ...
                final_R2, final_adjR2, final_p, pred_str);
        end
    end
end

%% 6) Summary Comparison Table
fprintf('\n============================================\n');
fprintf('=== SUMMARY: PE vs PCIst ===\n');
fprintf('============================================\n\n');

fprintf('%-8s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
    'Measure', 'Metric', '3M rho', '3M p', '6M rho', '6M p', '6M p_boot');
fprintf('%-8s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
    '--------', '----------', '----------', '----------', '----------', '----------', '----------');

for m = 1:length(measures)
    measure = measures{m};
    for k = 1:length(metrics)
        metric = metrics{k};
        r3 = corr_results.([measure '_' metric '_3M']);
        r6 = corr_results.([measure '_' metric '_6M']);
        pb = bootstrap_results.([measure '_' metric '_6M']);
        
        fprintf('%-8s %-12s %-12.3f %-12.4f %-12.3f %-12.4f %-12.4f\n', ...
            measure, strrep(metric,'_','-'), r3.rho, r3.p, r6.rho, r6.p, pb);
    end
end

%% 7) Save Results
save([directory, RESULTS_Folder, 'PE_PCIst_Statistics.mat'], ...
    'corr_results', 'bootstrap_results', 'reg_results', 'PE', 'PCIst', ...
    'GOSE_3M', 'GOSE_6M', 'subjects', 'raw_labels');

fid = fopen([directory, RESULTS_Folder, 'PE_PCIst_Summary.txt'], 'w');
fprintf(fid, 'PE and PCIst Statistical Analysis\n');
fprintf(fid, '=================================\n\n');
fprintf(fid, 'Subjects: %d (6-month: %d)\n\n', n_subjects, length(subs_6M));

fprintf(fid, 'CORRELATIONS:\n');
for m = 1:length(measures)
    measure = measures{m};
    for k = 1:length(metrics)
        metric = metrics{k};
        r3 = corr_results.([measure '_' metric '_3M']);
        r6 = corr_results.([measure '_' metric '_6M']);
        fprintf(fid, '  %s %s: 3M rho=%.3f p=%.4f | 6M rho=%.3f p=%.4f\n', ...
            measure, strrep(metric,'_','-'), r3.rho, r3.p, r6.rho, r6.p);
    end
end

fprintf(fid, '\nBOOTSTRAP:\n');
for m = 1:length(measures)
    measure = measures{m};
    for k = 1:length(metrics)
        metric = metrics{k};
        fprintf(fid, '  %s %s: 3M p_boot = %.4f | 6M p_boot = %.4f\n', ...
            measure, strrep(metric,'_','-'), ...
            bootstrap_results.([measure '_' metric '_3M']), ...
            bootstrap_results.([measure '_' metric '_6M']));
    end
end

fprintf(fid, '\nREGRESSION:\n');
for m = 1:length(measures)
    measure = measures{m};
    for k = 1:length(metrics)
        metric = metrics{k};
        for t = 1:2
            if t == 1, tp = '3M'; else, tp = '6M'; end
            r = reg_results.([measure '_' metric '_' tp]);
            fprintf(fid, '  %s %s %s: R2=%.3f, adjR2=%.3f, p=%.4f, retained: %s\n', ...
                measure, strrep(metric,'_','-'), tp, r.R2, r.adjR2, r.p, ...
                strjoin(r.predictors, ', '));
        end
    end
end

fclose(fid);

warning('on', 'stats:pca:ColRankDefX');

fprintf('\n=== STATISTICAL ANALYSIS COMPLETE ===\n');
fprintf('Results saved in: %s%s\n', directory, RESULTS_Folder);
fprintf('Figures saved in: %s%s\n', directory, FIGURES_Folder);