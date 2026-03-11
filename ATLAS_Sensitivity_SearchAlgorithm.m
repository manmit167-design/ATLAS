%% ATLAS: Sensitivity Analysis + Two-Stage Best-Fit Search
% Authors:
% 1) Dr. Manmit K Singh- Post Doctoral Fellow, 
% Department of Earth and Environmental Geosciences,
% University of Dayton, 300 College Park, Dayton, OH 45469
% 2) Dr. Umesh K Haritasya, 
% Director, Sustainability Program and Graduate Certificate
% Professor, Earth and Environmental Geosciences
% University of Dayton, 300 College Park, Dayton, OH 45469

%  Part A: Monte Carlo Sensitivity Analysis (LHS, Pearson r, S4.1-S4.3)
%  Part B: Two-Stage Hierarchical Bayesian Inversion (S5.1-S5.3)
%           Stage 1 – Sobol quasi-random global search (N=400)
%           Stage 2 – Nelder-Mead local refinement with Tikhonov regularisation


clc; clear; close all;

% ============================================================================
% Set random seed for reproducible sensitivity analysis
% ============================================================================
rng(42, 'twister');
fprintf('\n=== REPRODUCIBLE MODE ENABLED ===\n');
fprintf('Random seed set to: 42 (Mersenne Twister algorithm)\n');
fprintf('Sensitivity analysis weights will be identical across runs.\n\n');

% --------------------------- Globals for AOP_Rrs ---------------------------
global C_ph lambda S_CDOM type_Rrs_below zB ...
       type_case_water fA g_dd g_dsr g_dsa ...
       f_dd f_ds alpha view view_w sun sun_w ...
       rho_L P RH Hoz WV Rrs_obs

%% =========================================================================
%% PART A — SENSITIVITY ANALYSIS
%  Source: published sensitivity-only script.  Nothing changed.
%% =========================================================================

%% Step 1: Load Wavelength Data
[filename, pathname] = uigetfile('*.txt', 'Select the Wavelength Text File');
if isequal(filename, 0)
    error('User canceled file selection.');
end

filepath = fullfile(pathname, filename);
raw = load(filepath);

% Expecting: col1 = wavelength [nm], col2:end = satellite/reference spectra
lambda_full   = raw(:, 1);
Satellite_Ref = raw(:, 2:end);

% Sort if needed
if any(diff(lambda_full) < 0)
    warning('Wavelengths are not sorted. Sorting them.');
    [lambda_full, sort_idx] = sort(lambda_full);
    Satellite_Ref = Satellite_Ref(sort_idx, :);
end

% Show wavelengths with serial index
fprintf('\nAvailable Wavelengths:\n');
for ii = 1:length(lambda_full)
    fprintf('%3d : %.2f nm\n', ii, lambda_full(ii));
end

% Ask user for index range
start_idx = input('\nEnter the START serial number: ');
end_idx   = input('Enter the END serial number: ');

if start_idx < 1 || end_idx > length(lambda_full) || start_idx >= end_idx
    error('Invalid index range.');
end

% Selected working wavelength window
lambda = lambda_full(start_idx:end_idx);
n_wl   = numel(lambda);

% Dummy for AOP_Rrs
Rrs_obs = zeros(size(lambda));

%% Step 2: Monte Carlo Controls
N_MC = 100;
n_params = 5;   % [sun, view, CDOM, C_X, g_size]

% Fixed angular parameters
sun_mu_base  =  84.013;%EMIT% 20; %PRISMA 64;%Hyperion    % degrees
view_mu_base =  eps;%EMIT %20;    %PRISMA7;%Hyperion          % degrees

% Angular uncertainties (normal distribution)
sun_sigma_abs  = 5;   % degrees absolute
view_sigma_abs = 2;   % degrees absolute

% Parameter bounds
CDOM_bounds   = [0.01, 50];
CDOM_mu_geo   = 5.0;
CDOM_CV       = 0.5;

CX_bounds     = [0.1, 300];
CX_mu_geo     = 50;
CX_CV         = 0.7;

gsize_bounds  = [0.1e-6, 50e-6];
gsize_mu_geo  = 3.25e-6;
gsize_CV      = 0.4;

fprintf('\n=== PARAMETER BOUNDS ===\n');
fprintf('CDOM:      [%.2f, %.2f] mg/m³\n', CDOM_bounds(1), CDOM_bounds(2));
fprintf('SPM:       [%.2f, %.2f] g/m³\n', CX_bounds(1), CX_bounds(2));
fprintf('Grain Size: [%.2f, %.2f] μm\n', gsize_bounds(1)*1e6, gsize_bounds(2)*1e6);

% Convert to log-normal parameters  (CV-based — published formula)
[CDOM_mu_log, CDOM_sigma_log]   = geo_to_lognormal(CDOM_mu_geo, CDOM_CV);
[CX_mu_log, CX_sigma_log]       = geo_to_lognormal(CX_mu_geo, CX_CV);
[gsize_mu_log, gsize_sigma_log] = geo_to_lognormal(gsize_mu_geo, gsize_CV);

%% Wavelength-dependent weighting
lambda_weights = calculate_wavelength_weights(lambda);

%% ---------- Sensitivity Analysis with Latin Hypercube Sampling ----------
fprintf('\n========== Running Sensitivity Analysis (REPRODUCIBLE MODE) ==========\n');

fprintf('Generating Latin Hypercube Samples with fixed seed (42)...\n');
LHS_samples = lhsdesign(N_MC, n_params);

% Allocate arrays
Rrs_MC = zeros(n_wl, N_MC);
params = zeros(N_MC, n_params);

% Transform LHS samples to parameter distributions
sun_vals   = norminv(LHS_samples(:,1), sun_mu_base, sun_sigma_abs);
view_vals  = norminv(LHS_samples(:,2), view_mu_base, view_sigma_abs);

CDOM_vals  = logninv(LHS_samples(:,3), CDOM_mu_log, CDOM_sigma_log);
CX_vals    = logninv(LHS_samples(:,4), CX_mu_log, CX_sigma_log);
gsize_vals = logninv(LHS_samples(:,5), gsize_mu_log, gsize_sigma_log);

% Apply physical bounds
CDOM_vals  = max(CDOM_bounds(1),  min(CDOM_bounds(2),  CDOM_vals));
CX_vals    = max(CX_bounds(1),    min(CX_bounds(2),    CX_vals));
gsize_vals = max(gsize_bounds(1), min(gsize_bounds(2), gsize_vals));

% Run forward model
fprintf('Running %d Monte Carlo simulations...\n', N_MC);
for i = 1:N_MC
    if mod(i, 10) == 0
        fprintf('  Completed %d/%d\n', i, N_MC);
    end

    sun    = sun_vals(i);
    view   = view_vals(i);
    C_CDOM = CDOM_vals(i);
    C_X    = CX_vals(i);
    g_size = gsize_vals(i);

    S_CDOM = 0.0088 + 0.0092 * exp(-1 * C_CDOM);
    [view_w, sun_w, rho_L] = Snell_law(view, sun);

    C_ph = eps;
    type_Rrs_below = 0; zB = 4.0; type_case_water = 2;
    fA = [0,0,1,0,0,0]; g_dd = 0.05; g_dsr = 0; g_dsa = 0;
    f_dd = 1; f_ds = 1;
    alpha = 1.317; P = 1013.25; RH = 0.60; Hoz = 0.3; WV = 2.5;

    Fit = [C_CDOM, C_X, g_size];
    [Rrs0, ~] = AOP_Rrs_Mie(Fit);
    Rrs_MC(:, i) = real(Rrs0(:,1));
    params(i, :) = [sun, view, C_CDOM, C_X, g_size];
end

% Compute statistics
Rrs_mean_diag = mean(Rrs_MC, 2);
Rrs_std_diag  = std(Rrs_MC, 0, 2);
Rrs_CV_diag   = Rrs_std_diag ./ (abs(Rrs_mean_diag) + eps);

% Weighted correlation analysis
corr_matrix = zeros(n_params, n_wl);
for p = 1:n_params
    for w = 1:n_wl
        corr_matrix(p, w) = corr(params(:,p), Rrs_MC(w,:)', 'Type','Pearson', 'Rows','complete');
    end
end

% Weighted parameter importance
weighted_corr = abs(corr_matrix) .* repmat(lambda_weights', n_params, 1);
mean_weighted_corr = sum(weighted_corr, 2) / sum(lambda_weights);
weights = mean_weighted_corr / sum(mean_weighted_corr);

param_names = {'Sun Angle [deg]', 'View Angle [deg]', 'CDOM [mg/m^3]', 'C_X [g/m^3]', 'Grain Size [m]'};

fprintf('\n--- Weighted Parameter Importance (REPRODUCIBLE) ---\n');
for i = 1:n_params
    fprintf('%-20s : %.6f\n', param_names{i}, weights(i));
end

%% ---------- Visualization: 3-Panel Figure (Single Row) ----------
figure('Color','w','Position',[50 50 1800 500]);

% Panel 1: Rrs Mean ± 1σ
subplot(1,3,1);
hold on;
fill([lambda; flipud(lambda)], ...
     [Rrs_mean_diag - Rrs_std_diag; flipud(Rrs_mean_diag + Rrs_std_diag)], ...
     [0.8, 0.8, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(lambda, Rrs_mean_diag, 'b-', 'LineWidth', 2);
xlabel('Wavelength (nm)', 'FontSize', 12, 'FontWeight','bold');
ylabel('R_{rs} (sr^{-1})', 'FontSize', 12, 'FontWeight','bold');
set(gca, 'FontSize', 10, 'Box', 'off', 'LineWidth', 1.5);
ax1 = gca;
ax1.XAxis.LineWidth = 1.5;
ax1.YAxis.LineWidth = 1.5;
xlim([400 2500]);

% Panel 2: Parameter Sensitivity (Correlation)
subplot(1,3,2);
colors = lines(n_params);
hold on;
for p = 1:n_params
    plot(lambda, corr_matrix(p,:), 'LineWidth', 2, 'Color', colors(p,:));
end
xlabel('Wavelength (nm)', 'FontSize', 12, 'FontWeight','bold');
ylabel('Pearson r', 'FontSize', 12, 'FontWeight','bold');
legend(param_names, 'FontSize', 8, 'Location', 'best');
set(gca, 'FontSize', 10, 'Box', 'off', 'LineWidth', 1.5);
ax2 = gca;
ax2.XAxis.LineWidth = 1.5;
ax2.YAxis.LineWidth = 1.5;
xlim([400 2500]);
ylim([-1 1]);

% Panel 3: Weighted Parameter Importance
subplot(1,3,3);
bar(weights, 'FaceColor', [0.2 0.4 0.6]);
set(gca, 'XTickLabel', {'Sun','View','CDOM','C_X','g_size'}, ...
    'XTickLabelRotation', 45, 'FontSize', 10, 'Box', 'off', 'LineWidth', 1.5);
ylabel('Weighted Importance', 'FontSize', 12, 'FontWeight','bold');
ylim([0, max(weights)*1.15]);
ax3 = gca;
ax3.XAxis.LineWidth = 1.5;
ax3.YAxis.LineWidth = 1.5;

%% Save sensitivity results
save('Sensitivity_Analysis_Results.mat', 'lambda', 'corr_matrix', 'weights', ...
    'Rrs_mean_diag', 'Rrs_std_diag', 'Rrs_CV_diag', ...
    'params', 'Rrs_MC', 'CDOM_bounds', 'CX_bounds', 'gsize_bounds', ...
    'CDOM_CV', 'CX_CV', 'gsize_CV', 'lambda_weights', ...
    'LHS_samples', 'sun_vals', 'view_vals', 'CDOM_vals', 'CX_vals', 'gsize_vals');

fprintf('\n=============================================================\n');
fprintf('Sensitivity Analysis Complete!\n');
fprintf('Results saved to: Sensitivity_Analysis_Results.mat\n');
fprintf('=============================================================\n');
fprintf('\nREPRODUCIBILITY SUMMARY:\n');
fprintf('✓ Sensitivity analysis weights: FIXED (seed=42)\n');
fprintf('✓ Correlation matrix: REPRODUCIBLE\n');
fprintf('✓ Parameter importance ranking: CONSTANT\n');
fprintf('\nRun this code again → IDENTICAL results!\n');
fprintf('=============================================================\n');

%% =========================================================================
%% USER PERMISSION GATE
%% =========================================================================
fprintf('\n');
fprintf('══════════════════════════════════════════════════════════════\n');
fprintf('  Sensitivity analysis figure is ready.\n');
fprintf('  Proceed to two-stage best-fit search? (S5.1-S5.3)\n');
fprintf('  This will fit each reference spectrum in your data file.\n');
fprintf('  Stage 1: Sobol global search  (N_coarse = 400 evaluations)\n');
fprintf('  Stage 2: Nelder-Mead local refinement with Tikhonov prior\n');
fprintf('══════════════════════════════════════════════════════════════\n');
user_ans = input('  Enter  y  to proceed, any other key to exit: ', 's');
user_ans = strtrim(lower(user_ans));

if ~strcmp(user_ans, 'y')
    fprintf('\nExiting after sensitivity analysis. No best-fit search run.\n');
    return
end

%% =========================================================================
%% PART B — TWO-STAGE BEST-FIT SEARCH (S5.1-S5.3)
%  Source: original combined script, with one correction noted above.
%% =========================================================================
fprintf('\n========== Two-Stage Best-Fit Search (S5.1-S5.3) ==========\n');
fprintf('Using varied random sampling for thorough exploration.\n\n');

% Reset RNG to shuffle so optimization uses varied sampling (not fixed seed)
rng('shuffle');

% Reference spectra matched to selected wavelength window
Satellite_Ref_sel = Satellite_Ref(start_idx:end_idx, :);
n_refs = size(Satellite_Ref_sel, 2);

% Prior parameters (S5.2 — from glacial lake literature)
prior_means = [log(5.0); log(50); log(3.25e-6)];   % [CDOM; C_X; g_size]
prior_stds  = [0.50; 0.60; 0.40];                    % log-normal std deviations

% Optimisation settings (S5.3)
N_coarse   = 400;     % Sobol evaluations (Stage 1, S24)
lambda_reg = 0.01;    % Tikhonov regularisation strength (S5.2)
opts = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-7, ...
                'MaxIter', 500, 'MaxFunEvals', 2000);

% Output arrays
best_mu_all     = zeros(n_refs, 3);
best_rmse_all   = inf(n_refs, 1);
best_params_log = zeros(n_refs, 3);

% Pre-generate Sobol sequence once (same sequence for all reference spectra,
% ensures fair comparison across spectra — Stage 1 geometry is fixed)
sobol_seq     = sobolset(3);
sobol_samples = net(sobol_seq, N_coarse);   % N_coarse × 3, in [0,1]^3

% Total work units for progress bar
%   Stage 1: N_coarse evaluations per spectrum (trackable exactly)
%   Stage 2: ~500 iterations per spectrum (approximate; shown as indeterminate)
total_stage1_work = n_refs * N_coarse;

wb = waitbar(0, sprintf('Starting best-fit search (0/%d spectra)', n_refs), ...
    'Name', 'Two-Stage Best-Fit Search — S5.1-S5.3', ...
    'CreateCancelBtn', 'setappdata(gcbf,''cancelled'',true)');
setappdata(wb, 'cancelled', false);
t_start = tic;

for iRef = 1:n_refs

    % Check for user cancellation
    if getappdata(wb, 'cancelled')
        fprintf('\nUser cancelled best-fit search after spectrum %d.\n', iRef-1);
        break
    end

    target = Satellite_Ref_sel(:, iRef);
    fprintf('\n--- Fitting Reference Spectrum %d / %d ---\n', iRef, n_refs);

    % ── STAGE 1: Sobol Global Search (S5.1, Eq. S24-S26) ─────────────────
    fprintf('  Stage 1: Sobol global search (%d evaluations) ...\n', N_coarse);

    best_rmse   = inf;
    best_params = [NaN, NaN, NaN];
    valid_count = 0;

    for k = 1:N_coarse

        % Update progress bar every 20 samples
        if mod(k, 20) == 0 || k == 1
            work_done   = (iRef-1)*N_coarse + k;
            frac        = work_done / total_stage1_work;
            t_elapsed   = toc(t_start);
            if frac > 0
                t_remain = t_elapsed/frac * (1-frac);
                eta_str  = sprintf('ETA: %dm %02ds', floor(t_remain/60), mod(round(t_remain),60));
            else
                eta_str = 'ETA: calculating...';
            end
            waitbar(frac, wb, sprintf('Spectrum %d/%d  |  Stage 1: %d/%d  |  %s', ...
                iRef, n_refs, k, N_coarse, eta_str));
        end

        % Transform Sobol point to log-parameter space (Eq. S25)
        mu_CDOM_log  = log(CDOM_bounds(1)) + sobol_samples(k,1) * ...
                       (log(CDOM_bounds(2)) - log(CDOM_bounds(1)));
        mu_CX_log    = log(CX_bounds(1))   + sobol_samples(k,2) * ...
                       (log(CX_bounds(2))   - log(CX_bounds(1)));
        mu_gsize_log = log(gsize_bounds(1)) + sobol_samples(k,3) * ...
                       (log(gsize_bounds(2)) - log(gsize_bounds(1)));

        mu_CDOM  = exp(mu_CDOM_log);
        mu_CX    = exp(mu_CX_log);
        mu_gsize = exp(mu_gsize_log);

        try
            Rrs_try = run_mc_Rrs_mean_optimized(lambda, N_MC, ...
                sun_mu_base, view_mu_base, sun_sigma_abs, view_sigma_abs, ...
                mu_CDOM, mu_CX, mu_gsize, ...
                CDOM_sigma_log, CX_sigma_log, gsize_sigma_log);

            if any(~isfinite(Rrs_try)) || any(isnan(Rrs_try)) || any(Rrs_try <= 0)
                continue
            end

            % Weighted RMSE (Eq. S26)
            rmse_w = sqrt(sum(lambda_weights .* (target - Rrs_try).^2) / sum(lambda_weights));

            if ~isfinite(rmse_w) || isnan(rmse_w); continue; end

            valid_count = valid_count + 1;
            if rmse_w < best_rmse
                best_rmse   = rmse_w;
                best_params = [mu_CDOM, mu_CX, mu_gsize];
            end

        catch
            continue
        end
    end

    fprintf('  Stage 1 complete: %d valid / %d evaluated  |  best RMSE = %.4g\n', ...
        valid_count, N_coarse, best_rmse);

    % Fallback if no valid sample found
    if any(isnan(best_params)) || isinf(best_rmse)
        warning('Stage 1 found no valid sample. Falling back to prior means.');
        best_params = [exp(prior_means(1)), exp(prior_means(2)), exp(prior_means(3))];
        best_rmse   = inf;
    end

    fprintf('  Stage 1 best guess: CDOM=%.3f mg/m³, SPM=%.3f g/m³, g_size=%.3f μm\n', ...
        best_params(1), best_params(2), best_params(3)*1e6);

    % ── STAGE 2: Nelder-Mead Local Refinement (S5.1, Eq. S27-S29) ────────
    fprintf('  Stage 2: Nelder-Mead refinement (λ_reg=%.2f) ...\n', lambda_reg);
    waitbar(((iRef-1)*N_coarse + N_coarse)/total_stage1_work, wb, ...
        sprintf('Spectrum %d/%d  |  Stage 2: Nelder-Mead refining...', iRef, n_refs));

    obj_fun = @(x_log) objective_function_bayesian_safe(x_log, lambda, N_MC, target, ...
        sun_mu_base, view_mu_base, sun_sigma_abs, view_sigma_abs, ...
        CDOM_sigma_log, CX_sigma_log, gsize_sigma_log, ...
        CDOM_bounds, CX_bounds, gsize_bounds, lambda_weights, ...
        prior_means, prior_stds, lambda_reg);

    x0_log = log(best_params(:));
    if any(~isfinite(x0_log)) || any(isnan(x0_log))
        warning('Stage 1 initial guess invalid in log-space. Using prior means.');
        x0_log = prior_means;
    end

    try
        [x_opt_log, fval, exitflag] = fminsearch(obj_fun, x0_log, opts);
        if exitflag <= 0
            warning('Nelder-Mead did not converge (exitflag=%d). Keeping Stage 1 result.', exitflag);
            x_opt_log = x0_log;
            fval      = best_rmse;
        end
    catch ME
        warning('Nelder-Mead failed: %s. Keeping Stage 1 result.', ME.message);
        x_opt_log = x0_log;
        fval      = best_rmse;
    end

    best_params_opt = exp(x_opt_log);
    best_params_opt(1) = max(CDOM_bounds(1),  min(CDOM_bounds(2),  best_params_opt(1)));
    best_params_opt(2) = max(CX_bounds(1),    min(CX_bounds(2),    best_params_opt(2)));
    best_params_opt(3) = max(gsize_bounds(1), min(gsize_bounds(2), best_params_opt(3)));

    if any(isnan(best_params_opt)) || any(~isfinite(best_params_opt))
        warning('Stage 2 returned invalid values. Keeping Stage 1 result.');
        best_params_opt = best_params;
        fval            = best_rmse;
    end

    % Accept Stage 2 result only if it is strictly better than Stage 1
    if isfinite(fval) && fval < best_rmse
        best_mu_all(iRef,:)      = best_params_opt(:)';
        best_rmse_all(iRef)      = fval;
        best_params_log(iRef,:)  = log(best_params_opt(:))';
    else
        best_mu_all(iRef,:)      = best_params(:)';
        best_rmse_all(iRef)      = best_rmse;
        best_params_log(iRef,:)  = log(best_params(:))';
    end

    % ── Per-spectrum results ──────────────────────────────────────────────
    fprintf('\n  *** FINAL BEST-FIT  (Spectrum %d) ***\n', iRef);
    fprintf('    CDOM       = %.3f mg/m³\n',  best_mu_all(iRef,1));
    fprintf('    C_X (SPM)  = %.3f g/m³\n',   best_mu_all(iRef,2));
    fprintf('    Grain Size = %.3f μm  (%.3e m)\n', ...
            best_mu_all(iRef,3)*1e6, best_mu_all(iRef,3));
    fprintf('    Weighted RMSE = %.4g\n', best_rmse_all(iRef));

    % 95% Confidence Intervals (±2 log-normal std dev, S5.3)
    CDOM_CI  = exp([log(best_mu_all(iRef,1)) - 2*CDOM_sigma_log, ...
                    log(best_mu_all(iRef,1)) + 2*CDOM_sigma_log]);
    CX_CI    = exp([log(best_mu_all(iRef,2)) - 2*CX_sigma_log, ...
                    log(best_mu_all(iRef,2)) + 2*CX_sigma_log]);
    gsize_CI = exp([log(best_mu_all(iRef,3)) - 2*gsize_sigma_log, ...
                    log(best_mu_all(iRef,3)) + 2*gsize_sigma_log]);

    fprintf('    95%% CI:\n');
    fprintf('      CDOM:      [%.3f, %.3f] mg/m³\n', CDOM_CI(1),      CDOM_CI(2));
    fprintf('      C_X:       [%.3f, %.3f] g/m³\n',  CX_CI(1),        CX_CI(2));
    fprintf('      g_size:    [%.3f, %.3f] μm\n',    gsize_CI(1)*1e6, gsize_CI(2)*1e6);

    if best_mu_all(iRef,3) < 1e-6
        warning('Grain size < 1 μm — may indicate insufficient SWIR coverage or model limits.');
    elseif best_mu_all(iRef,3) > 20e-6
        warning('Grain size > 20 μm — may indicate non-glacial sediment contribution.');
    else
        fprintf('    ✓ Grain size within expected glacial flour range (1–20 μm)\n');
    end
end

% Close waitbar
if ishandle(wb); delete(wb); end
t_total = toc(t_start);
fprintf('\nTotal elapsed time: %dm %02ds\n', floor(t_total/60), mod(round(t_total),60));

%% ── Final Summary Table ──────────────────────────────────────────────────
fprintf('\n========== FINAL SUMMARY ==========\n');
fprintf('%-6s %-14s %-14s %-16s %-12s\n', ...
    'Ref#', 'CDOM[mg/m³]', 'SPM[g/m³]', 'g_size[μm]', 'RMSE');
fprintf('%-6s %-14s %-14s %-16s %-12s\n', ...
    '----', '------------', '----------', '-----------', '----');
for i = 1:n_refs
    fprintf('%-6d %-14.3f %-14.3f %-16.3f %-12.4g\n', ...
        i, best_mu_all(i,1), best_mu_all(i,2), best_mu_all(i,3)*1e6, best_rmse_all(i));
end

if n_refs > 1
    fprintf('\nStatistical Summary (all reference spectra):\n');
    fprintf('  CDOM:       Mean=%.3f  Std=%.3f  Range=[%.3f, %.3f] mg/m³\n', ...
        mean(best_mu_all(:,1)), std(best_mu_all(:,1)), ...
        min(best_mu_all(:,1)),  max(best_mu_all(:,1)));
    fprintf('  SPM:        Mean=%.3f  Std=%.3f  Range=[%.3f, %.3f] g/m³\n', ...
        mean(best_mu_all(:,2)), std(best_mu_all(:,2)), ...
        min(best_mu_all(:,2)),  max(best_mu_all(:,2)));
    fprintf('  Grain Size: Mean=%.3f  Std=%.3f  Range=[%.3f, %.3f] μm\n', ...
        mean(best_mu_all(:,3))*1e6, std(best_mu_all(:,3))*1e6, ...
        min(best_mu_all(:,3))*1e6,  max(best_mu_all(:,3))*1e6);
end

%% Save combined results
save('GLAM_BioLithRT_Sensitivity_BestFit_Results.mat', ...
    'lambda', 'corr_matrix', 'weights', ...
    'Rrs_mean_diag', 'Rrs_std_diag', 'Rrs_CV_diag', ...
    'params', 'Rrs_MC', 'lambda_weights', ...
    'CDOM_bounds', 'CX_bounds', 'gsize_bounds', ...
    'CDOM_CV', 'CX_CV', 'gsize_CV', ...
    'LHS_samples', 'sun_vals', 'view_vals', 'CDOM_vals', 'CX_vals', 'gsize_vals', ...
    'best_mu_all', 'best_rmse_all', 'best_params_log', ...
    'prior_means', 'prior_stds', 'N_coarse', 'lambda_reg');

fprintf('\nResults saved: GLAM_BioLithRT_Sensitivity_BestFit_Results.mat\n');
fprintf('=============================================================\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function [mu_log, sigma_log] = geo_to_lognormal(mu_geo, CV)
    % Convert geometric mean and CV to log-normal parameters.
    % Published formula — do not change.
    sigma_log = sqrt(log(1 + CV^2));
    mu_log    = log(mu_geo);
end

function weights = calculate_wavelength_weights(lambda)
    % Uniform weighting across VIS/NIR/SWIR — published formula, do not change.
    weights = ones(size(lambda));
    vis_idx  = lambda >= 400  & lambda < 700;   weights(vis_idx)  = 1.0;
    nir_idx  = lambda >= 700  & lambda < 1000;  weights(nir_idx)  = 1.0;
    swir_idx = lambda >= 1000 & lambda <= 2500; weights(swir_idx) = 1.0;
    weights  = weights / mean(weights);
end

function Rrs_mean = run_mc_Rrs_mean_optimized(lambda, N_MC, sun_mu, view_mu, ...
                                               sun_sigma, view_sigma, ...
                                               CDOM_mu, CX_mu, gsize_mu, ...
                                               CDOM_sigma_log, CX_sigma_log, gsize_sigma_log)
    % Forward model used by the optimiser (Stage 1 evaluation + Stage 2 objective).
    %
    % S_CDOM CORRECTION: original code used 0.001 + 0.008*exp(-C_CDOM).
    % Corrected here to 0.0088 + 0.0092*exp(-C_CDOM) to match the published
    % sensitivity forward model.  Using two different S_CDOM formulas in
    % sensitivity vs optimisation is an internal inconsistency that would
    % cause the optimizer to search a landscape computed by a different model
    % than the one the sensitivity analysis characterised.

    global C_ph S_CDOM type_Rrs_below zB type_case_water fA g_dd g_dsr g_dsa ...
           f_dd f_ds alpha view view_w sun sun_w rho_L P RH Hoz WV

    n_wl_local = numel(lambda);
    Rrs_stack  = zeros(n_wl_local, N_MC);

    CDOM_mu_log  = log(CDOM_mu);
    CX_mu_log    = log(CX_mu);
    gsize_mu_log = log(gsize_mu);

    C_ph            = eps;
    type_Rrs_below  = 0;
    zB              = 4.0;
    type_case_water = 2;
    fA              = [0,0,1,0,0,0];
    g_dd = 0.05; g_dsr = 0; g_dsa = 0;
    f_dd = 1;    f_ds  = 1;
    alpha = 1.317; P = 1013.25; RH = 0.60; Hoz = 0.3; WV = 2.5;

    for k = 1:N_MC
        sun  = normrnd(sun_mu,  sun_sigma);
        view = normrnd(view_mu, view_sigma);

        C_CDOM = lognrnd(CDOM_mu_log,  CDOM_sigma_log);
        C_X    = lognrnd(CX_mu_log,    CX_sigma_log);
        g_size = lognrnd(gsize_mu_log, gsize_sigma_log);

        % Corrected S_CDOM — consistent with published sensitivity model
        S_CDOM = 0.0088 + 0.0092 * exp(-1 * C_CDOM);

        [view_w, sun_w, rho_L] = Snell_law(view, sun);

        Fit = [C_CDOM, C_X, g_size];
        [Rrs0, ~]        = AOP_Rrs_Mie(Fit);
        Rrs_stack(:, k)  = real(Rrs0(:,1));
    end

    Rrs_mean = mean(Rrs_stack, 2);
end

function cost = objective_function_bayesian_safe(x_log, lambda, N_MC, target, ...
                                                  sun_mu, view_mu, sun_sigma, view_sigma, ...
                                                  CDOM_sigma_log, CX_sigma_log, gsize_sigma_log, ...
                                                  CDOM_bounds, CX_bounds, gsize_bounds, ...
                                                  lambda_weights, prior_means, prior_stds, lambda_reg)
    % Bayesian objective function with Tikhonov regularisation (S5.2).
    % cost = sqrt(chi2_data) + lambda_reg * sqrt(chi2_prior)
    % All NaN/Inf guards retained from original.

    x_log = x_log(:);

    if any(~isfinite(x_log)) || any(isnan(x_log))
        cost = 1e10; return
    end

    params_lin = exp(x_log);

    if params_lin(1) < CDOM_bounds(1) || params_lin(1) > CDOM_bounds(2) || ...
       params_lin(2) < CX_bounds(1)   || params_lin(2) > CX_bounds(2)   || ...
       params_lin(3) < gsize_bounds(1)|| params_lin(3) > gsize_bounds(2)
        cost = 1e10; return
    end

    try
        Rrs_pred = run_mc_Rrs_mean_optimized(lambda, N_MC, ...
            sun_mu, view_mu, sun_sigma, view_sigma, ...
            params_lin(1), params_lin(2), params_lin(3), ...
            CDOM_sigma_log, CX_sigma_log, gsize_sigma_log);

        if any(~isfinite(Rrs_pred)) || any(isnan(Rrs_pred)) || any(Rrs_pred <= 0)
            cost = 1e10; return
        end
    catch
        cost = 1e10; return
    end

    residuals = target - Rrs_pred;
    if any(~isfinite(residuals)) || any(isnan(residuals))
        cost = 1e10; return
    end

    chi2_data = sum(lambda_weights .* residuals.^2) / sum(lambda_weights);
    if ~isfinite(chi2_data) || isnan(chi2_data) || chi2_data < 0
        cost = 1e10; return
    end

    prior_means = prior_means(:);
    prior_stds  = prior_stds(:);
    chi2_prior  = sum(((x_log - prior_means) ./ prior_stds).^2);
    if ~isfinite(chi2_prior) || isnan(chi2_prior) || chi2_prior < 0
        cost = 1e10; return
    end

    cost = sqrt(chi2_data) + lambda_reg * sqrt(chi2_prior);

    % Boundary penalty — softly discourages solutions near hard bounds
    penalty = 0;
    if params_lin(1) < CDOM_bounds(1)*1.05 || params_lin(1) > CDOM_bounds(2)*0.95
        penalty = penalty + 0.2*cost;
    end
    if params_lin(2) < CX_bounds(1)*1.05   || params_lin(2) > CX_bounds(2)*0.95
        penalty = penalty + 0.2*cost;
    end
    if params_lin(3) < gsize_bounds(1)*1.1 || params_lin(3) > gsize_bounds(2)*0.9
        penalty = penalty + 0.3*cost;
    end
    cost = cost + penalty;

    if ~isfinite(cost) || isnan(cost) || cost < 0
        cost = 1e10;
    end

    cost = cost(1);
end