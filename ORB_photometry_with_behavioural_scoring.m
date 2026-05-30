%% ========================================================================
%  FIBER PHOTOMETRY ANALYSIS COUPLED WITH BEHAVIOURAL SCORING
%  ORB kisspeptin neuronal activity

% =========================================================================

clear; close all; clc;

%% ========================================================================
%  SECTION 1: CONFIGURATION
% =========================================================================
PRE_SEC  = 5;           % seconds before event to extract
POST_SEC = 3;           % seconds after event to extract
BL_START = -5;          % baseline window start (for z-scoring)
BL_END   = -1;          % baseline window end
PERI_START = 0;         % peri-onset analysis window start
PERI_END   = 3;         % peri-onset analysis window end
N_PERM   = 1000;        % number of permutation iterations
RESP_THRESH = 1.0;      % z-score threshold for "responsive" bout

%% ========================================================================
%  SECTION 2: DATA LOADING
%  
%  dFF is pre-computed as: (GCaMP_465 - Fitted_405_IRLS) / Fitted_405_IRLS
%  The 405nm isosbestic was fitted using iteratively reweighted least 
%  squares (IRLS) to remove motion artifacts and photobleaching.
%  Sampling rate: ~12.05 Hz (dt ~ 0.083 s) across all recordings.
% =========================================================================

% Define conditions and their files
conditions = struct();

% --- Condition 1: Object exploration ---
conditions(1).name = 'Object exploration';
conditions(1).short = 'object';
conditions(1).animals = {'4764','5122','5123','5124','5313','5314'};
conditions(1).color = [0.902 0.318 0.000];  % orange
conditions(1).file_pattern = '%s_object_doric_IRLS_Output_withbehavior.xlsx';

% --- Condition 2: Male partner ---
conditions(2).name = 'Male partner';
conditions(2).short = 'male';
conditions(2).animals = {'4764','5122','5123','5124','5314'};
conditions(2).color = [0.416 0.106 0.604];  % purple
conditions(2).file_pattern = '%s_male_doric_IRLS_Output_withbehavior.xlsx';

% --- Condition 3: Female partner (sexual behavior) ---
conditions(3).name = 'Female partner';
conditions(3).short = 'female';
conditions(3).animals = {'4764','5122','5123','5124','5314'};
conditions(3).color = [0.130 0.400 0.675];  % blue
conditions(3).file_pattern = '%s_female_doric_IRLS_Output_withbehavior.xlsx';

% --- Load all data ---
data = struct();
for ci = 1:length(conditions)
    cond = conditions(ci);
    for ai = 1:length(cond.animals)
        animal = cond.animals{ai};
            
            filename = sprintf(cond.file_pattern, animal);
            T = readtable(filename);
            
            time_vec = T.Time_s;
            dff_vec  = T.dFF;
            
            % behavior_type column has START/STOP strings mixed with NaN
            % In MATLAB, NaN cells may read as missing/empty
            beh_col = T.behavior_type;
            if iscell(beh_col)
                idx_start = strcmp(beh_col, 'START');
                idx_stop  = strcmp(beh_col, 'STOP');
            else
                % If read as categorical
                idx_start = beh_col == "START";
                idx_stop  = beh_col == "STOP";
            end
            starts = T.Time(idx_start);
            stops  = T.Time(idx_stop);
      %  end
        
        data(ci).animal{ai}  = animal;
        data(ci).time{ai}    = time_vec;
        data(ci).dff{ai}     = dff_vec;
        data(ci).starts{ai}  = starts;
        data(ci).stops{ai}   = stops;
        data(ci).dt(ai)      = median(diff(time_vec));
    end
end

fprintf('Data loaded successfully.\n');
for ci = 1:length(conditions)
    total_bouts = 0;
    for ai = 1:length(conditions(ci).animals)
        total_bouts = total_bouts + length(data(ci).starts{ai});
    end
    fprintf('  %s: n=%d animals, %d bouts\n', ...
        conditions(ci).name, length(conditions(ci).animals), total_bouts);
end

%% ========================================================================
%  SECTION 3: CORE ANALYSIS FUNCTIONS
%
%  These two functions form the backbone of the entire analysis.
%  Every figure and test uses them.
% =========================================================================

% --- Function 1: Extract peri-event traces ---
% For each event time, cut a window from the continuous signal.
% Inputs:
%   time   - [Nx1] time vector from recording
%   dff    - [Nx1] dFF signal
%   events - [Mx1] event times (behavior onsets or offsets)
%   pre    - seconds before event (positive number)
%   post   - seconds after event (positive number)
% Outputs:
%   t_axis - [1 x n_samples] time axis relative to event
%   traces - [n_valid_events x n_samples] matrix of extracted traces
%
% Algorithm:
%   1. Compute dt = median(diff(time))
%   2. n_pre = round(pre/dt), n_post = round(post/dt)
%   3. For each event, find closest index in time vector
%   4. Extract dff(idx-n_pre : idx+n_post-1)
%   5. Discard if window falls outside recording boundaries

% --- Function 2: Z-score traces to pre-event baseline ---
% Each trace is independently normalized using its own baseline.
% Inputs:
%   t_axis - [1 x T] time axis
%   traces - [N x T] matrix of traces
%   bl_start, bl_end - baseline window boundaries (e.g., -5, -1)
% Output:
%   z_traces - [N x T] z-scored traces
%
% Algorithm:
%   1. Find columns where bl_start <= t_axis < bl_end
%   2. For each row: mu = mean of baseline samples, sd = std of baseline
%   3. z = (trace - mu) / sd
%   4. If sd < 1e-8, set sd = 1e-8 (avoid division by zero)
%
% NOTE: Because z-scoring uses the Pre window as baseline, the Pre window
% values will always center near zero by construction. The statistical
% question is whether Peri and Post windows deviate from this baseline.

%% ========================================================================
%  SECTION 4: RUN ANALYSIS FOR ALL CONDITIONS
% =========================================================================

results = struct();

for ci = 1:length(conditions)
    n_animals = length(conditions(ci).animals);
    
    all_z_onset  = [];   % pooled z-scored onset traces (all bouts)
    all_z_offset = [];   % pooled z-scored offset traces
    animal_peri  = zeros(1, n_animals);  % per-animal peri-onset mean z
    animal_traces_onset = cell(1, n_animals);  % per-animal z traces
    animal_traces_offset = cell(1, n_animals);  % per-animal z traces
    
    for ai = 1:n_animals
        time_vec = data(ci).time{ai};
        dff_vec  = data(ci).dff{ai};
        starts   = data(ci).starts{ai};
        stops    = data(ci).stops{ai};
        dt       = data(ci).dt(ai);
        
        % --- Extract onset-aligned traces ---
        [t_ax, traces_on] = extract_peri_event(time_vec, dff_vec, starts, PRE_SEC, POST_SEC);
        z_on = zscore_traces(t_ax, traces_on, BL_START, BL_END);
        
        all_z_onset = [all_z_onset; z_on];
        animal_traces_onset{ai} = z_on;
        
        % Compute peri-onset mean for this animal
        peri_mask = t_ax >= PERI_START & t_ax < PERI_END;
        animal_peri(ai) = mean(mean(z_on(:, peri_mask), 2));
        
        % --- Extract offset-aligned traces ---
        [~, traces_off] = extract_peri_event(time_vec, dff_vec, stops, PRE_SEC, POST_SEC);
        z_off = zscore_traces(t_ax, traces_off, BL_START, BL_END);
        all_z_offset = [all_z_offset; z_off];
        animal_traces_offset{ai} = z_off;
    end
    
    % Store results
    results(ci).t_axis = t_ax;
    results(ci).all_z_onset  = all_z_onset;
    results(ci).all_z_offset = all_z_offset;
    results(ci).animal_peri  = animal_peri;
    results(ci).animal_traces_onset = animal_traces_onset;
    results(ci).animal_traces_offset = animal_traces_offset;
    results(ci).n_bouts = size(all_z_onset, 1);
    results(ci).n_animals = n_animals;
end

%% ========================================================================
%  SECTION 5: STATISTICAL TESTS (run for each condition)
% =========================================================================

for ci = 1:length(conditions)
    fprintf('\n=== %s (n=%d animals, %d bouts) ===\n', ...
        conditions(ci).name, results(ci).n_animals, results(ci).n_bouts);
    
    ap = results(ci).animal_peri;
    t_ax = results(ci).t_axis;
    
    % --- TEST 1: One-sample t-test (peri-onset z vs zero) ---
    % Tests whether the group mean peri-onset z differs from zero.
    % Each animal contributes one value (its mean peri-onset z across bouts).
    [~, p_ttest, ~, stats_t] = ttest(ap);
    fprintf('  One-sample t-test: t=%.3f, p=%.4f\n', stats_t.tstat, p_ttest);
    
    % --- TEST 2: Wilcoxon signed-rank test (one-sided, greater) ---
    % Non-parametric test appropriate for small n.
    % Tests whether the distribution of animal peri-onset z values is
    % shifted above zero.
    % NOTE: With n=5, minimum possible one-sided p = 1/2^5 = 0.0312
    %        With n=4, minimum p = 0.0625 (significance at alpha=0.05 
    %        is mathematically impossible)
    p_wilcox = signrank(ap, 0, 'tail', 'right');
    fprintf('  Wilcoxon signed-rank (one-sided): p=%.4f\n', p_wilcox);
    
    % --- TEST 3: Permutation test ---
    % Distribution-free test of whether neural activity is time-locked
    % to behavior events.
    %
    % Algorithm:
    %   1. Compute real group mean peri-onset z (average across animals)
    %   2. For each of N_PERM iterations:
    %      a. For each animal, generate random event times uniformly
    %         distributed within the valid recording window
    %         [time(1)+PRE+1, time(end)-POST-1]
    %         Generate same number of fake events as real bouts
    %      b. Extract peri-event traces at these random times
    %      c. Z-score using same baseline window
    %      d. Compute mean peri-onset z for this animal
    %      e. Average across animals -> one shuffled group mean
    %   3. p-value = proportion of shuffled means >= real mean
    
    n_animals = results(ci).n_animals;
    real_means = zeros(1, n_animals);
    shuffled_group = zeros(N_PERM, 1);
    
    for ai = 1:n_animals
        time_vec = data(ci).time{ai};
        dff_vec  = data(ci).dff{ai};
        starts   = data(ci).starts{ai};
        
        % Real peri-onset z
        [t_s, tr] = extract_peri_event(time_vec, dff_vec, starts, PRE_SEC, POST_SEC);
        z = zscore_traces(t_s, tr, BL_START, BL_END);
        pm = t_s >= PERI_START & t_s < PERI_END;
        real_means(ai) = mean(mean(z(:, pm), 2));
        
        % Shuffled
        t_min = time_vec(1) + PRE_SEC + 1;
        t_max = time_vec(end) - POST_SEC - 1;
        n_events = length(starts);
        
        shuf_means = zeros(N_PERM, 1);
        for p = 1:N_PERM
            fake_starts = t_min + (t_max - t_min) * rand(n_events, 1);
            [ts, trs] = extract_peri_event(time_vec, dff_vec, fake_starts, PRE_SEC, POST_SEC);
            zs = zscore_traces(ts, trs, BL_START, BL_END);
            shuf_means(p) = mean(mean(zs(:, pm), 2));
        end
        % Store per-animal shuffled means (will average across animals)
        if ai == 1
            all_shuf = zeros(N_PERM, n_animals);
        end
        all_shuf(:, ai) = shuf_means;
    end
    
    real_group = mean(real_means);
    shuffled_group = mean(all_shuf, 2);  % average across animals per iteration
    p_perm = mean(shuffled_group >= real_group);
    
    fprintf('  Permutation test (%d iter): real z=%.4f, p=%.4f\n', ...
        N_PERM, real_group, p_perm);
    
    % --- TEST 4: Paired t-tests between time windows ---
    % (Reported for Figure 1 panels G and H)
    % For each animal: compute mean z in Pre, Peri, Post windows
    animal_pre  = zeros(1, n_animals);
    animal_post = zeros(1, n_animals);
    for ai = 1:n_animals
        z_traces = results(ci).animal_traces_onset{ai};
        pre_mask  = t_ax >= -5 & t_ax < -1;
        peri_mask = t_ax >= 0  & t_ax < 3;
        post_mask = t_ax >= 3  & t_ax < 8;
        animal_pre(ai)  = mean(mean(z_traces(:, pre_mask), 2));
        results(ci).animal_peri(ai) = mean(mean(z_traces(:, peri_mask), 2));
        animal_post(ai) = mean(mean(z_traces(:, post_mask), 2));
    end
    [~, p_pre_peri] = ttest(animal_pre, results(ci).animal_peri);
    [~, p_pre_post] = ttest(animal_pre, animal_post);
    fprintf('  Paired t-test Pre vs Peri: p=%.4f\n', p_pre_peri);
    fprintf('  Paired t-test Pre vs Post: p=%.4f\n', p_pre_post);
    
    % --- TEST 5: AUC analysis ---
    % Area under curve of mean z-scored trace per animal
    % Pre-onset AUC: trapz over [-5, -1]s
    % Post-onset AUC: trapz over [0, 5]s
    pre_auc_mask  = t_ax >= -5 & t_ax < -1;
    post_auc_mask = t_ax >= 0  & t_ax < 3;
    pre_auc  = zeros(1, n_animals);
    post_auc = zeros(1, n_animals);
    for ai = 1:n_animals
        mean_trace = mean(results(ci).animal_traces_onset{ai}, 1);
        pre_auc(ai)  = trapz(t_ax(pre_auc_mask), mean_trace(pre_auc_mask));
        post_auc(ai) = trapz(t_ax(post_auc_mask), mean_trace(post_auc_mask));
    end
    [~, p_auc] = ttest(pre_auc, post_auc);
    fprintf('  AUC onset Pre vs Post: p=%.4f\n', p_auc);
    
    % --- TEST 6: Bout responsiveness ---
    % Percentage of bouts with mean peri-onset z > threshold
    n_resp = 0; n_total = 0;
    resp_per_animal = zeros(1, n_animals);
    for ai = 1:n_animals
        z_traces = results(ci).animal_traces_onset{ai};
        peri_vals = mean(z_traces(:, peri_mask), 2);
        resp_per_animal(ai) = 100 * sum(peri_vals > RESP_THRESH) / length(peri_vals);
        n_resp  = n_resp + sum(peri_vals > RESP_THRESH);
        n_total = n_total + length(peri_vals);
    end
    fprintf('  Responsive bouts: %.1f%% (range: %.0f-%.0f%%)\n', ...
        mean(resp_per_animal), min(resp_per_animal), max(resp_per_animal));
    
    % Store everything
    results(ci).p_ttest  = p_ttest;
    results(ci).p_wilcox = p_wilcox;
    results(ci).p_perm   = p_perm;
    results(ci).real_group = real_group;
    results(ci).shuffled_group = shuffled_group;
    results(ci).resp_pct = mean(resp_per_animal);
    results(ci).resp_per_animal = resp_per_animal;
    results(ci).pre_auc  = pre_auc;
    results(ci).post_auc = post_auc;
    results(ci).p_auc    = p_auc;
    results(ci).animal_pre  = animal_pre;
    results(ci).animal_post = animal_post;
end

%% ========================================================================
%  SECTION 6: FIGURE 1 - Female Partner PETH (8-panel)
%
%  This is the primary finding figure. 8 panels (A-H):
%    A: Grand average onset-aligned trace (all bouts pooled)
%    B: Grand average offset-aligned trace
%    C: Single-bout heatmap (onset), sorted by peak time
%    D: Single-bout heatmap (offset)
%    E: Per-animal average onset traces
%    F: Per-animal average offset traces
%    G: Onset window quantification (Pre/Peri/Post bars)
%    H: Offset window quantification
% =========================================================================

% % ci = 2;  % Female partner
% % r = results(ci);
% % t = r.t_axis;
% % 
% % figure('Position', [50 50 1400 1600], 'Color', 'w');
% % 
% % % --- Panel A: Grand average onset ---
% % subplot(4,2,1);
% % mu = mean(r.all_z_onset, 1);
% % se = std(r.all_z_onset, 0, 1) / sqrt(size(r.all_z_onset, 1));
% % fill([t fliplr(t)], [mu-se fliplr(mu+se)], conditions(ci).color, ...
% %     'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
% % plot(t, mu, 'Color', conditions(ci).color, 'LineWidth', 2);
% % xline(0, '--k', 'LineWidth', 1.5);
% % yline(0, ':k', 'Alpha', 0.5);
% % xlabel('Time from behavior onset (s)');
% % ylabel('z-scored \DeltaF/F');
% % title('A. Onset-aligned (all bouts pooled)', 'FontWeight', 'bold');
% % xlim([-PRE_SEC POST_SEC]);
% % 
% % % --- Panel B: Grand average offset ---
% % subplot(4,2,2);
% % mu_off = mean(r.all_z_offset, 1);
% % se_off = std(r.all_z_offset, 0, 1) / sqrt(size(r.all_z_offset, 1));
% % fill([t fliplr(t)], [mu_off-se_off fliplr(mu_off+se_off)], [0.694 0.094 0.169], ...
% %     'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
% % plot(t, mu_off, 'Color', [0.694 0.094 0.169], 'LineWidth', 2);
% % xline(0, '--k', 'LineWidth', 1.5);
% % yline(0, ':k', 'Alpha', 0.5);
% % xlabel('Time from behavior offset (s)');
% % ylabel('z-scored \DeltaF/F');
% % title('B. Offset-aligned (all bouts pooled)', 'FontWeight', 'bold');
% % xlim([-PRE_SEC POST_SEC]);
% % 
% % % --- Panel C: Heatmap onset ---
% % subplot(4,2,3);
% % post_cols = t >= 0;
% % [~, peak_idx] = max(r.all_z_onset(:, post_cols), [], 2);
% % [~, sort_order] = sort(peak_idx);
% % imagesc(t, 1:size(r.all_z_onset,1), r.all_z_onset(sort_order, :));
% % caxis([-3 3]); colormap(redblue_colormap()); colorbar;
% % hold on; xline(0, '--w', 'LineWidth', 1.5);
% % xlabel('Time from onset (s)');
% % ylabel('Bout (sorted)');
% % title('C. Single-bout heatmap (onset)', 'FontWeight', 'bold');
% % 
% % % --- Panel D: Heatmap offset ---
% % subplot(4,2,4);
% % [~, peak_idx_off] = max(r.all_z_offset(:, post_cols), [], 2);
% % [~, sort_off] = sort(peak_idx_off);
% % imagesc(t, 1:size(r.all_z_offset,1), r.all_z_offset(sort_off, :));
% % caxis([-3 3]); colorbar;
% % hold on; xline(0, '--w', 'LineWidth', 1.5);
% % xlabel('Time from offset (s)');
% % ylabel('Bout (sorted)');
% % title('D. Single-bout heatmap (offset)', 'FontWeight', 'bold');
% % 
% % % --- Panel E: Per-animal onset ---
% % subplot(4,2,5);
% % colors_a = lines(length(conditions(ci).animals));
% % for ai = 1:length(conditions(ci).animals)
% %     m = mean(r.animal_traces_onset{ai}, 1);
% %     plot(t, m, 'Color', colors_a(ai,:), 'LineWidth', 1.5); hold on;
% % end
% % xline(0, '--k'); yline(0, ':k');
% % xlabel('Time from onset (s)');
% % ylabel('z-scored \DeltaF/F');
% % title('E. Per-animal averages (onset)', 'FontWeight', 'bold');
% % legend(conditions(ci).animals, 'Location', 'northeast');
% % xlim([-PRE_SEC POST_SEC]);
% % 
% % % --- Panel F: Per-animal offset ---
% % subplot(4,2,6);
% % for ai = 1:length(conditions(ci).animals)
% %     time_vec = data(ci).time{ai};
% %     dff_vec  = data(ci).dff{ai};
% %     stops    = data(ci).stops{ai};
% %     [~, tr_off] = extract_peri_event(time_vec, dff_vec, stops, PRE_SEC, POST_SEC);
% %     z_off = zscore_traces(t, tr_off, BL_START, BL_END);
% %     plot(t, mean(z_off, 1), 'Color', colors_a(ai,:), 'LineWidth', 1.5); hold on;
% % end
% % xline(0, '--k'); yline(0, ':k');
% % xlabel('Time from offset (s)');
% % ylabel('z-scored \DeltaF/F');
% % title('F. Per-animal averages (offset)', 'FontWeight', 'bold');
% % legend(conditions(ci).animals, 'Location', 'northeast');
% % xlim([-PRE_SEC POST_SEC]);
% % 
% % % --- Panel G: Onset bars ---
% % subplot(4,2,7);
% % bar_data = [mean(r.animal_pre), mean(r.animal_peri), mean(r.animal_post)];
% % bar_err  = [std(r.animal_pre)/sqrt(length(r.animal_pre)), ...
% %             std(r.animal_peri)/sqrt(length(r.animal_peri)), ...
% %             std(r.animal_post)/sqrt(length(r.animal_post))];
% % b = bar(1:3, bar_data); hold on;
% % errorbar(1:3, bar_data, bar_err, 'k.', 'LineWidth', 1.5);
% % % Individual dots
% % for ai = 1:length(r.animal_pre)
% %     plot(1 + 0.1*randn, r.animal_pre(ai), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
% %     plot(2 + 0.1*randn, r.animal_peri(ai), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
% %     plot(3 + 0.1*randn, r.animal_post(ai), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
% % end
% % set(gca, 'XTickLabel', {'Pre (-5:-1s)', 'Peri (0:3s)', 'Post (3:8s)'});
% % ylabel('z-scored \DeltaF/F');
% % title('G. Onset quantification', 'FontWeight', 'bold');
% % 
% % % --- Panel H: Offset bars (same logic, different windows) ---
% % subplot(4,2,8);
% % % Offset windows: Pre -5:-1, Peri -1:1, Post 1:8
% % off_pre_mask  = t >= -5 & t < -1;
% % off_peri_mask = t >= -1 & t < 1;
% % off_post_mask = t >= 1  & t < 8;
% % aoff_pre = zeros(1, results(ci).n_animals);
% % aoff_peri = zeros(1, results(ci).n_animals);
% % aoff_post = zeros(1, results(ci).n_animals);
% % for ai = 1:results(ci).n_animals
% %     time_vec = data(ci).time{ai};
% %     dff_vec  = data(ci).dff{ai};
% %     stops    = data(ci).stops{ai};
% %     [~, tr_off] = extract_peri_event(time_vec, dff_vec, stops, PRE_SEC, POST_SEC);
% %     z_off = zscore_traces(t, tr_off, BL_START, BL_END);
% %     aoff_pre(ai)  = mean(mean(z_off(:, off_pre_mask), 2));
% %     aoff_peri(ai) = mean(mean(z_off(:, off_peri_mask), 2));
% %     aoff_post(ai) = mean(mean(z_off(:, off_post_mask), 2));
% % end
% % bar_off = [mean(aoff_pre), mean(aoff_peri), mean(aoff_post)];
% % bar_off_err = [std(aoff_pre)/sqrt(length(aoff_pre)), ...
% %                std(aoff_peri)/sqrt(length(aoff_peri)), ...
% %                std(aoff_post)/sqrt(length(aoff_post))];
% % bar(1:3, bar_off); hold on;
% % errorbar(1:3, bar_off, bar_off_err, 'k.', 'LineWidth', 1.5);
% % for ai = 1:length(aoff_pre)
% %     plot(1+0.1*randn, aoff_pre(ai), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
% %     plot(2+0.1*randn, aoff_peri(ai), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
% %     plot(3+0.1*randn, aoff_post(ai), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
% % end
% % set(gca, 'XTickLabel', {'Pre (-5:-1s)', 'Peri (-1:1s)', 'Post (1:8s)'});
% % ylabel('z-scored \DeltaF/F');
% % title('H. Offset quantification', 'FontWeight', 'bold');
% % 
% % sgtitle(sprintf('Figure 1: Neural Activity During Female Partner Interaction\nn=%d animals, %d bouts', ...
% %     results(ci).n_animals, results(ci).n_bouts), 'FontWeight', 'bold', 'FontSize', 14);
% % 
% % saveas(gcf, 'Figure1_FemalePartner_PETH.png');

%% ========================================================================
%  SECTION 7: FIGURE 2 - Three-Condition Comparison
%
%  Three rows:
%    Row 1 (A-D): Grand average PETH per condition
%    Row 2 (E-H): Single-bout heatmaps per condition
%    Row 3: (I) Cross-condition bar plot, (J) Summary stats table
% =========================================================================

figure('Position', [50 50 1600 1200], 'Color', 'w');

for ci = 1:3
    r = results(ci);
    t = r.t_axis;
    
    % Row 1: PETH traces
    subplot(3, 3, ci);
    mu = mean(r.all_z_onset, 1);
    se = std(r.all_z_onset, 0, 1) / sqrt(size(r.all_z_onset, 1));
    fill([t fliplr(t)], [mu-se fliplr(mu+se)], conditions(ci).color, ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
    plot(t, mu, 'Color', conditions(ci).color, 'LineWidth', 2);
    xline(0, '--k'); yline(0, ':k');
    xlim([-PRE_SEC POST_SEC]);
    xticks(-5:1:3);
    ylim([-0.3 0.4]);
    yticks(-0.3:0.1:0.4);
    xlabel('Time (s)'); 
    if ci == 1; ylabel('z-scored \DeltaF/F'); end
    pm = t >= 0 & t < 3;
    title(sprintf('%s\n(n=%d, %d bouts)\nz=%.3f', conditions(ci).name, ...
        r.n_animals, r.n_bouts, mean(r.animal_peri)), ...
        'FontWeight', 'bold', 'Color', conditions(ci).color);
    
    % Row 2: Heatmaps
%     subplot(3, 3, 3+ci);
%     post_cols = t >= 0;
%     [~, pk] = max(r.all_z_onset(:, post_cols), [], 2);
%     [~, so] = sort(pk);
%     imagesc(t, 1:r.n_bouts, r.all_z_onset(so, :));
%      n = 256;
%     cmap = [linspace(0,1,n/2)' linspace(0,1,n/2)' ones(n/2,1); ...
%         ones(n/2,1) linspace(1,0,n/2)' linspace(1,0,n/2)'];
% colormap(cmap)
%     caxis([-3 3]); hold on; xline(0, '--w');
%     xlim([-5 3]);   % restrict x-axis
%     hold on; xline(0, '--w');
%     xlabel('Time (s)');
%     if ci == 1; ylabel('Bout (sorted)'); end
%     title(sprintf('n = %d bouts', r.n_bouts), 'Color', [0.5 0.5 0.5]);
%     cb = colorbar;
%     cb.Label.String = 'z-scored \DeltaF/F';

 
% Row 3: Permutations
    subplot(3, 3, 6+ci);
    % subplot(3, 3, (ci-1)*3 + 3);
    histogram(r.shuffled_group, 35, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', [0.6 0.6 0.6]);
    hold on; xline(r.real_group, '--r', 'LineWidth', 2);
    sig_str = '';
    if r.p_perm < 0.05; sig_str = ' *'; end
    if ci == 1; title(sprintf('Permutation p=%.3f%s', r.p_perm, sig_str), 'FontWeight', 'bold');
    else; title(sprintf('p=%.3f%s', r.p_perm, sig_str)); end

end

% % % Row 3 left: Bar plot
% % subplot(3, 4, [9 10]);
% % means_all = arrayfun(@(r) mean(r.animal_peri), results);
% % sems_all  = arrayfun(@(r) std(r.animal_peri)/sqrt(r.n_animals), results);
% % colors_bar = vertcat(conditions.color);
% % for ci = 1:3
% %     bar(ci, means_all(ci), 'FaceColor', colors_bar(ci,:), 'EdgeColor', 'k'); hold on;
% % end
% % errorbar(1:3, means_all, sems_all, 'k.', 'LineWidth', 1.5);
% % % Animal dots
% % for ci = 1:3
% %     ap = results(ci).animal_peri;
% %     scatter(ci + 0.1*randn(size(ap)), ap, 25, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
% % end
% % set(gca, 'XTick', 1:3, 'XTickLabel', {conditions.short});
% % ylabel('Post-onset z (0-3s)');
% % title('Post-onset activity across conditions', 'FontWeight', 'bold');
% % yline(0, ':k');

% % % Row 3 right: Stats table
% % subplot(3, 3, [11 12]); axis off;
% % col_headers = {'Condition', 'n', 'Bouts', 'Peri z', 'Perm p', 'Wilcox p', 'Resp %'};
% % table_data = cell(4, 7);
% % for ci = 1:3
% %     table_data{ci,1} = conditions(ci).short;
% %     table_data{ci,2} = num2str(results(ci).n_animals);
% %     table_data{ci,3} = num2str(results(ci).n_bouts);
% %     table_data{ci,4} = sprintf('%.3f', mean(results(ci).animal_peri));
% %     table_data{ci,5} = sprintf('%.3f', results(ci).p_perm);
% %     table_data{ci,6} = sprintf('%.3f', results(ci).p_wilcox);
% %     table_data{ci,7} = sprintf('%.1f', results(ci).resp_pct);
% % end
% % % Display as text
% % for row = 0:3
% %     for col = 1:7
% %         if row == 0
% %             txt = col_headers{col};
% %             fw = 'bold';
% %         else
% %             txt = table_data{row, col};
% %             fw = 'normal';
% %         end
% %         text(0.02 + (col-1)*0.14, 0.9 - row*0.18, txt, ...
% %             'FontSize', 9, 'FontWeight', fw, 'Units', 'normalized');
% %     end
% % end
% % title('Statistical comparison', 'FontWeight', 'bold');

sgtitle('Figure 2: Behavioral Selectivity Across Three Paradigms', ...
    'FontWeight', 'bold', 'FontSize', 14);

saveas(gcf, 'Figure2_ThreeCondition_Comparison.png');

% %% ========================================================================
% %  SECTION 8: SUPPLEMENTAL FIGURE S1 - Raw Traces
% %
% %  For each condition, plot raw dFF traces for all animals with
% %  behavior bouts overlaid as colored shading.
% % =========================================================================
% 
% figure('Position', [50 50 1800 2000], 'Color', 'w');
% panel = 0;
% total_panels = sum(cellfun(@length, {conditions.animals}));
% 
% for ci = 1:3
%     for ai = 1:length(conditions(ci).animals)
%         panel = panel + 1;
%         subplot(total_panels, 1, panel);
% 
%         plot(data(ci).time{ai}, data(ci).dff{ai}, 'Color', [0.2 0.2 0.2], 'LineWidth', 0.2);
%         hold on;
% 
%         % Overlay behavior bouts
%         starts = data(ci).starts{ai};
%         stops  = data(ci).stops{ai};
%         n_bouts = min(length(starts), length(stops));
%         for bi = 1:n_bouts
%             patch([starts(bi) stops(bi) stops(bi) starts(bi)], ...
%                   [min(ylim) min(ylim) max(ylim) max(ylim)], ...
%                   conditions(ci).color, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
%         end
% 
%         ylabel('\DeltaF/F', 'FontSize', 7);
%         title(sprintf('%s - Mouse %s (%d bouts)', conditions(ci).name, ...
%             conditions(ci).animals{ai}, length(starts)), ...
%             'FontSize', 8, 'FontWeight', 'bold', 'Color', conditions(ci).color);
%         set(gca, 'FontSize', 7);
%     end
% end
% xlabel('Time (s)');
% sgtitle('Supplemental Figure S1: Raw Traces', 'FontWeight', 'bold');
% saveas(gcf, 'FigureS1_RawTraces.png');

%% ========================================================================
%  SECTION 9: SUPPLEMENTAL FIGURE S2 - Detailed Stats (4x4 grid)
%
%  Rows = conditions, Columns = analyses:
%    Col 1: Per-animal onset traces
%    Col 2: AUC paired plot
%    Col 3: Permutation test histogram
%    Col 4: Bout responsiveness per animal
% =========================================================================

figure('Position', [50 50 1600 1200], 'Color', 'w');

for ci = 1:3
    r = results(ci);
    t = r.t_axis;
    colors_a = lines(r.n_animals);
    
    % Column 1: Per-animal traces
    subplot(4, 4, (ci-1)*3 + 1);
    for ai = 1:r.n_animals
        plot(t, mean(r.animal_traces_onset{ai}, 1), ...
            'Color', colors_a(ai,:), 'LineWidth', 1.2); hold on;
    end
    xline(0, '--k'); yline(0, ':k'); xlim([-PRE_SEC POST_SEC]);
    if ci == 1; title('Per-animal onset', 'FontWeight', 'bold'); end
    ylabel(conditions(ci).short, 'FontWeight', 'bold', 'Color', conditions(ci).color);
    
    % Column 2: AUC
    subplot(4, 4, (ci-1)*3 + 2);
    for ai = 1:r.n_animals
        plot([1 2], [r.pre_auc(ai) r.post_auc(ai)], 'o-', ...
            'Color', colors_a(ai,:), 'MarkerSize', 6, 'LineWidth', 1.2); hold on;
    end
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Pre','Post'});
    if ci == 1; title(sprintf('AUC (p=%.3f)', r.p_auc), 'FontWeight', 'bold');
    else; title(sprintf('p=%.3f', r.p_auc)); end
    
    % % Column 3: Permutation histogram
    % subplot(3, 3, (ci-1)*3 + 3);
    % histogram(r.shuffled_group, 35, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', [0.6 0.6 0.6]);
    % hold on; xline(r.real_group, '--r', 'LineWidth', 2);
    % sig_str = '';
    % if r.p_perm < 0.05; sig_str = ' *'; end
    % if ci == 1; title(sprintf('Permutation p=%.3f%s', r.p_perm, sig_str), 'FontWeight', 'bold');
    % else; title(sprintf('p=%.3f%s', r.p_perm, sig_str)); end
    
    % Column 4: Responsiveness
    subplot(4, 4, (ci-1)*3 + 4);
    bar(r.resp_per_animal, 'FaceColor', conditions(ci).color, 'EdgeColor', 'k');
    hold on; yline(mean(r.resp_per_animal), '--r');
    set(gca, 'XTickLabel', conditions(ci).animals, 'FontSize', 7);
    if ci == 1; title(sprintf('Resp. %.0f%%', r.resp_pct), 'FontWeight', 'bold');
    else; title(sprintf('%.0f%%', r.resp_pct)); end
end

sgtitle('Supplemental Figure S2: Detailed Statistics', 'FontWeight', 'bold');
saveas(gcf, 'FigureS2_DetailedStats.png');

fprintf('\n=== ALL FIGURES SAVED ===\n');

%% ========================================================================
%  HELPER FUNCTIONS (place these at the end of the script or in separate files)
% =========================================================================

function [t_axis, traces] = extract_peri_event(time, dff, events, pre_sec, post_sec)
    % EXTRACT_PERI_EVENT Extract peri-event traces from continuous signal
    %
    %   [t_axis, traces] = extract_peri_event(time, dff, events, pre_sec, post_sec)
    %
    %   time   - [Nx1] continuous time vector
    %   dff    - [Nx1] continuous dFF signal
    %   events - [Mx1] event times
    %   pre_sec, post_sec - window in seconds
    %
    %   Returns:
    %   t_axis - [1 x n_samples] time relative to event
    %   traces - [n_valid x n_samples] extracted traces
    
    dt = median(diff(time));
    n_pre  = round(pre_sec / dt);
    n_post = round(post_sec / dt);
    n_total = n_pre + n_post;
    t_axis = linspace(-pre_sec, post_sec, n_total);
    
    traces = [];
    for i = 1:length(events)
        [~, idx] = min(abs(time - events(i)));
        si = idx - n_pre;
        ei = idx + n_post - 1;
        if si >= 1 && ei <= length(dff) && (ei - si + 1) == n_total
            traces = [traces; dff(si:ei)'];
        end
    end
end

function z_traces = zscore_traces(t_axis, traces, bl_start, bl_end)
    % ZSCORE_TRACES Z-score each trace to its own pre-event baseline
    %
    %   z_traces = zscore_traces(t_axis, traces, bl_start, bl_end)
    %
    %   Each row is independently normalized:
    %     z = (trace - mean(baseline)) / std(baseline)
    %   where baseline = samples in [bl_start, bl_end)
    
    bl_mask = t_axis >= bl_start & t_axis < bl_end;
    bl_mean = mean(traces(:, bl_mask), 2);
    bl_std  = std(traces(:, bl_mask), 0, 2);
    bl_std(bl_std < 1e-8) = 1e-8;  % prevent division by zero
    z_traces = (traces - bl_mean) ./ bl_std;
end

function cmap = redblue_colormap()
    % REDBLUE_COLORMAP Generate red-white-blue diverging colormap
    n = 256;
    r = [linspace(0.019,1,n/2) linspace(1,0.694,n/2)];
    g = [linspace(0.188,1,n/2) linspace(1,0.094,n/2)];
    b = [linspace(0.380,1,n/2) linspace(1,0.169,n/2)];
    cmap = [r' g' b'];
end
% %% %% ========================================================================
% %  SECTION 10: Excel tables
% %
% % =========================================================================
% % % Prepare table for Excel
% 
% % --- Prepare onset-aligned Excel table with mouse and bout info ---
% t_axis = results(1).t_axis;           % time relative to onset
% z_traces = results(1).all_z_onset;    % pooled z-scored onset traces
% n_bouts = size(z_traces,1);
% 
% % Create labels for each bout: "MouseID_BoutX"
% mouse_labels = {};
% for ai = 1:results(1).n_animals
%     n_bouts_animal = size(results(1).animal_traces_onset{ai},1);
%     for bi = 1:n_bouts_animal
%         mouse_labels{end+1,1} = sprintf('%s_Bout%d', conditions(1).animals{ai}, bi);
%     end
% end
% 
% % Check length matches
% if length(mouse_labels) ~= n_bouts
%     warning('Number of labels does not match number of rows in all_z_onset');
% end
% 
% % Create table with numeric column headers
% var_names = compose('%.2f', t_axis);   % just numbers as strings
% T_excel = array2table(z_traces, 'VariableNames', var_names);
% 
% % Add Mouse/Bout column at the beginning
% T_excel = addvars(T_excel, mouse_labels, 'Before', 1, 'NewVariableNames', {'Mouse_Bout'});
% 
% %% --- Append per-mouse mean rows (mean across bouts) ---
% 
% mean_rows = [];
% mean_labels = {};
% 
% for ai = 1:results(1).n_animals
%     % Mean across all bouts for this mouse
%     mouse_mean = mean(results(1).animal_traces_onset{ai}, 1);
% 
%     mean_rows = [mean_rows; mouse_mean];
%     mean_labels{end+1,1} = sprintf('%s_MEAN', conditions(1).animals{ai});
% end
% 
% % Convert to table with same column names
% T_mean = array2table(mean_rows, 'VariableNames', var_names);
% 
% % Add Mouse label column
% T_mean = addvars(T_mean, mean_labels, ...
%     'Before', 1, 'NewVariableNames', {'Mouse_Bout'});
% 
% % Append to existing per-bout table
% T_excel = [T_excel; T_mean];
% 
% % Write Excel file
% writetable(T_excel, 'Figure1A_onset_aligned_perBout.xlsx');

%% ========================================================================
%  SECTION 10: Excel tables (ONSET, all conditions)
% =========================================================================

for ci = 1:length(conditions)

    % --- Grab condition-specific data ---
    t_axis   = results(ci).t_axis;          % time relative to onset
    z_traces = results(ci).all_z_onset;     % pooled z-scored onset traces
    n_bouts  = size(z_traces, 1);

    % --- Create labels for each bout: MouseID_BoutX ---
    mouse_labels = {};
    for ai = 1:results(ci).n_animals
        n_bouts_animal = size(results(ci).animal_traces_onset{ai}, 1);
        for bi = 1:n_bouts_animal
            mouse_labels{end+1,1} = sprintf('%s_Bout%d', ...
                conditions(ci).animals{ai}, bi);
        end
    end

    % Safety check
    if length(mouse_labels) ~= n_bouts
        warning('%s: Number of labels does not match number of bouts', ...
            conditions(ci).name);
    end

    % --- Create table ---
    var_names = compose('%.2f', t_axis);
    T_excel = array2table(z_traces, 'VariableNames', var_names);

    % Add Mouse/Bout column
    T_excel = addvars(T_excel, mouse_labels, ...
        'Before', 1, 'NewVariableNames', {'Mouse_Bout'});

    %% --- Append per-mouse mean rows ---
    mean_rows   = [];
    mean_labels = {};

    for ai = 1:results(ci).n_animals
        mouse_mean = mean(results(ci).animal_traces_onset{ai}, 1);
        mean_rows = [mean_rows; mouse_mean];
        mean_labels{end+1,1} = sprintf('%s_MEAN', ...
            conditions(ci).animals{ai});
    end

    T_mean = array2table(mean_rows, 'VariableNames', var_names);
    T_mean = addvars(T_mean, mean_labels, ...
        'Before', 1, 'NewVariableNames', {'Mouse_Bout'});

    % Append
    T_excel = [T_excel; T_mean];

    % --- Write file ---
    outname = sprintf('Figure1_onset_%s_perBout.xlsx', ...
        conditions(ci).short);
    writetable(T_excel, outname);

    fprintf('Saved onset table: %s\n', outname);

end

% %% offset 

% --- Prepare onset-aligned Excel table with mouse and bout info ---
t_axis = results(1).t_axis;           % time relative to offset
z_traces = results(1).all_z_offset;    % pooled z-scored offset traces
n_bouts = size(z_traces,1);

% Create labels for each bout: "MouseID_BoutX"
mouse_labels = {};
for ai = 1:results(1).n_animals
    n_bouts_animal = size(results(1).animal_traces_offset{ai},1);
    for bi = 1:n_bouts_animal
        mouse_labels{end+1,1} = sprintf('%s_Bout%d', conditions(1).animals{ai}, bi);
    end
end

% Check length matches
if length(mouse_labels) ~= n_bouts
    warning('Number of labels does not match number of rows in all_z_offset');
end

% Create table with numeric column headers
var_names = compose('%.2f', t_axis);   % just numbers as strings
T_excel = array2table(z_traces, 'VariableNames', var_names);

% Add Mouse/Bout column at the beginning
T_excel = addvars(T_excel, mouse_labels, 'Before', 1, 'NewVariableNames', {'Mouse_Bout'});

%% --- Append per-mouse mean rows (mean across bouts) ---

mean_rows = [];
mean_labels = {};

for ai = 1:results(1).n_animals
    % Mean across all bouts for this mouse
    mouse_mean = mean(results(1).animal_traces_offset{ai}, 1);

    mean_rows = [mean_rows; mouse_mean];
    mean_labels{end+1,1} = sprintf('%s_MEAN', conditions(1).animals{ai});
end

% Convert to table with same column names
T_mean = array2table(mean_rows, 'VariableNames', var_names);

% Add Mouse label column
T_mean = addvars(T_mean, mean_labels, ...
    'Before', 1, 'NewVariableNames', {'Mouse_Bout'});

% Append to existing per-bout table
T_excel = [T_excel; T_mean];

% Write Excel file
writetable(T_excel, 'Figure1B_offset_aligned_perBout.xlsx');
