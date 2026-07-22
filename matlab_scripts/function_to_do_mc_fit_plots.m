function [h, stats, handles, pcovar, limits] = ...
    function_to_do_mc_fit_plots( ...
    x, y, z, symb, fit, xname, yname, zname, sname, uname, model, filter, ...
    pcovar, limits, opts, subspec)
%
% Manipulate and plot the data in *_model.pts.
%
% JADC, 8/2025
% Refactored April 5, 2026

if nargin < 16 || isempty(opts)
    opts = struct();
end
if ~isfield(opts,'plotMPP'), opts.plotMPP = false; end
if ~isfield(opts,'plotcov'), opts.plotcov = true;  end
if ~isfield(opts,'ploterr'), opts.ploterr = true;  end
if ~isfield(opts,'plotbst'), opts.plotbst = false; end
if ~isfield(opts,'printz'),  opts.printz  = false; end
if ~isfield(opts,'print'),   opts.print   = false; end
if ~isfield(opts,'maksub'),        opts.maksub        = true;  end
if ~isfield(opts,'sigma_level'),   opts.sigma_level   = 1;     end
if ~isfield(opts,'legend_coord_digits'), opts.legend_coord_digits = 4; end
if ~isfield(opts,'legend_sigma_digits'), opts.legend_sigma_digits = 2; end
if ~isfield(opts,'cap_fraction'),  opts.cap_fraction  = 0.10;  end
if ~isfield(opts,'target_rse'),    opts.target_rse    = 0.20;  end
if ~isfield(opts,'best_override'), opts.best_override = []; end
if ~isfield(opts,'target_axes'), opts.target_axes = []; end
if ~isfield(opts,'filled_symbols_for_nofit'), opts.filled_symbols_for_nofit = true; end
if ~isfield(opts,'imprecision_filter_info'), opts.imprecision_filter_info = struct('applied',false); end
if ~isfield(opts,'central_epsz_override'), opts.central_epsz_override = NaN; end
opts.filled_symbols_for_nofit = logical(opts.filled_symbols_for_nofit);

sc = ['b','r','b','r','b','r','b','r','b','r','b','r'];
mt = ['o','o','s','s','v','v','d','d','s','s','v','v'];
logz = false;
dlw = get(groot,'DefaultLineLineWidth');

bayes = double(zname == "ln(Bayes)") * 3;
pert_sym = 2 + bayes;
cent_sym = 3 + bayes;

idx_best = find(symb == 1 & ~isnan(x) & ~isnan(y) & ~isnan(z));
idx_pert = find(symb == pert_sym & ~isnan(x) & ~isnan(y) & ~isnan(z));
idx_cent = find(symb == cent_sym & ~isnan(x) & ~isnan(y) & ~isnan(z));
idx_minrng = find(symb == 9 & ~isnan(x) & ~isnan(y) & ~isnan(z));
idx_maxrng = find(symb == 10 & ~isnan(x) & ~isnan(y) & ~isnan(z));

if model == "Central"
    % If the original best central model has already been suppressed upstream
    % (e.g. because the user chose to remove models that do not fit all data),
    % choose the best remaining central model from the finite central trials.
    if isempty(idx_best) && ~isempty(idx_cent)
        if bayes == 3
            [~, k] = max(z(idx_cent));
        else
            [~, k] = min(z(idx_cent));
        end
        idx_best = idx_cent(k);
    end
end

if isempty(idx_best)
    error('No best central model could be identified in %s data.', char(model))
end

idx_best = idx_best(1);
idx_cent(idx_cent == idx_best) = [];

idx_trials = idx_pert;
if model == "Central"
    idx_trials = idx_cent;
end

% Plotting/statistical set for the current analysis.
if model == "Perturbation"
    idx_sigma = idx_pert;
else
    idx_sigma = [idx_best; idx_cent(:)];
end
idx_sigma = idx_sigma(:).';
idx_sigma = idx_sigma(~isnan(idx_sigma));

% Finite rows actually used for covariance/imprecision estimates.  Outlier
% filtering is applied at read time by setting discarded perturbation rows to
% NaN, so this count is after z-filtering.
data_xyz_sigma = [x(idx_sigma), y(idx_sigma), z(idx_sigma)];
valid_xyz_sigma = all(isfinite(data_xyz_sigma), 2);
idx_sigma_cov = idx_sigma(valid_xyz_sigma);

% Global limits excluding marker 0, 9, 10 and NaNs.
valid_lim = ~isnan(x) & ~isnan(y) & ~isnan(z) & ~ismember(symb, [0 9 10]);
if any(valid_lim)
    limits.minx = min(limits.minx, min(x(valid_lim)));
    limits.maxx = max(limits.maxx, max(x(valid_lim)));
    limits.miny = min(limits.miny, min(y(valid_lim)));
    limits.maxy = max(limits.maxy, max(y(valid_lim)));
    limits.minz = min(limits.minz, min(z(valid_lim)));
    limits.maxz = max(limits.maxz, max(z(valid_lim)));
end

% MPP coordinate from symbols 9 and 10.
if ~isempty(idx_minrng) && ~isempty(idx_maxrng)
    MPPX = (x(idx_minrng(1)) + x(idx_maxrng(1))) / 2;
    MPPY = (y(idx_minrng(1)) + y(idx_maxrng(1))) / 2;
    MPPZ = (z(idx_minrng(1)) + z(idx_maxrng(1))) / 2;
else
    MPPX = NaN; MPPY = NaN; MPPZ = NaN;
end

stats = struct();
stats.model = model;
stats.best_idx = idx_best;
stats.best_x = x(idx_best);
stats.best_y = y(idx_best);
stats.best_z = z(idx_best);
stats.zopt = z(idx_best);
stats.best_fit = fit(idx_best);
stats.n_sigma_points = numel(idx_sigma_cov);
stats.idx_sigma = idx_sigma;
stats.idx_sigma_cov = idx_sigma_cov;
stats.sigx = NaN;
stats.sigy = NaN;
stats.epsz = NaN;
stats.filter_half_width_z = NaN;
stats.sigxt = NaN;
stats.sigyt = NaN;
stats.pcovar = [];
stats.imprecision_filter_info = opts.imprecision_filter_info;
stats.aleatoric_equated_to_data = false;
stats.aleatoric_covariance_available = true;
stats.overall_idx = NaN;
stats.overall_x = NaN;
stats.overall_y = NaN;
stats.overall_z = NaN;
stats.overall_fit = NaN;

if ~isempty(opts.best_override)
    stats.best_x = opts.best_override.x;
    stats.best_y = opts.best_override.y;
    stats.best_z = opts.best_override.z;
    if isfield(opts.best_override,'fit') && ~isempty(opts.best_override.fit)
        stats.best_fit = opts.best_override.fit;
    end
end

handles = struct('best', gobjects(1,1), ...
    'trials_fit', gobjects(1,1), ...
    'trials_nofit', gobjects(1,1), ...
    'overall', gobjects(1,1), ...
    'mpp', gobjects(1,1), ...
    'cov', gobjects(1,1), ...
    'cov_total', gobjects(1,1), ...
    'data_cov', gobjects(1,1), ...
    'data', [], ...
    'imprecision', [], ...
    'total', []);

if ~isempty(opts.target_axes) && isgraphics(opts.target_axes,'axes')
    h(1) = opts.target_axes;
    axes(h(1)); %#ok<LAXES>
elseif opts.maksub
    h(1) = subplot(subspec{1}{:});
else
    figure(1)
    h(1) = gca;
end
hold(h(1),'on')

% Best central model marker.
if ~isempty(opts.best_override)
    handles.best = plot_best_model_coords(h(1), stats.best_x, stats.best_y, stats.best_z, stats.best_fit, dlw, opts.filled_symbols_for_nofit);
else
    handles.best = plot_best_model(h(1), x, y, z, fit, idx_best, dlw, opts.filled_symbols_for_nofit);
end

% Error bars and statistics.
if opts.ploterr
    if model == "Perturbation"
        [handles.data, stats] = draw_perturbation_errorbars(h(1), x, y, z, idx_sigma, stats, dlw, opts.sigma_level, opts.cap_fraction);
    else
        [handles.imprecision, handles.total, stats] = draw_central_errorbars( ...
            h(1), x, y, z, idx_sigma, idx_best, pcovar, filter, stats, dlw, opts.sigma_level, opts.cap_fraction, opts.central_epsz_override);
    end
else
    stats = compute_stats_only(x, y, z, idx_sigma, stats, model, filter, idx_best, pcovar);
end

if ~isempty(opts.best_override)
    stats.best_x = opts.best_override.x;
    stats.best_y = opts.best_override.y;
    stats.best_z = opts.best_override.z;
    if isfield(opts.best_override,'fit') && ~isempty(opts.best_override.fit)
        stats.best_fit = opts.best_override.fit;
    end
end

% Trial points.
if ~isempty(idx_trials)
    isymb = symb(idx_trials(1));
    [handles.trials_fit, handles.trials_nofit, h_dum] = ...
        plot_trials(h(1), x, y, z, fit, idx_trials, sc(isymb), mt(isymb), dlw, opts.filled_symbols_for_nofit);
else
    h_dum = plot(nan,nan,'LineStyle','none');
end

% Best overall model.  In a perturbation plot this is the best finite
% perturbation model, but it is only plotted/reported if it is better than
% the best central model shown in the same plot.  If no perturbation model
% improves on the central optimum, the best central model is already the
% overall best model and no separate "Best overall model" marker is drawn.
if opts.plotbst && model == "Perturbation" && ~isempty(idx_pert)
    zpert = z(idx_pert);
    if bayes == 3
        [~, kp] = max(zpert);
    else
        [~, kp] = min(zpert);
    end

    idx_overall_candidate = idx_pert(kp);
    if bayes == 3
        candidate_is_better = z(idx_overall_candidate) > stats.best_z;
    else
        candidate_is_better = z(idx_overall_candidate) < stats.best_z;
    end

    if candidate_is_better
        stats.overall_idx = idx_overall_candidate;
        stats.overall_x = x(stats.overall_idx);
        stats.overall_y = y(stats.overall_idx);
        stats.overall_z = z(stats.overall_idx);
        stats.overall_fit = fit(stats.overall_idx);
        handles.overall = plot_best_overall(h(1), x, y, z, fit, stats.overall_idx, dlw, opts.filled_symbols_for_nofit);
    end
end

format_axes(h(1), xname, yname, zname, logz)

% Covariance sample size gate and covariance ellipse(s).
axes(h(1))
hold(h(1),'on')
n_cov_points = numel(idx_sigma_cov);

n_req = ceil(1 + 1/(2*opts.target_rse^2));
handled_sample_warning = false;

if model == "Central" && filter && isfield(opts,'imprecision_filter_info') && ...
        isstruct(opts.imprecision_filter_info) && isfield(opts.imprecision_filter_info,'applied') && ...
        opts.imprecision_filter_info.applied

    info = opts.imprecision_filter_info;
    n_pre = info.pre_count;
    n_post = info.post_count;
    n_removed = info.removed_count;
    handled_sample_warning = true;

    if isfinite(n_pre) && n_pre < n_req
        fpre = ceil(n_req / max(n_pre,1));
        msg = sprintf(['Only %d central model(s) were available before imprecision filtering; ' ...
            'at least %d are required for the estimated standard deviation to be accurate to within %.0f%%.\n\n' ...
            'This is a pre-filtering sample-size problem rather than evidence for negligible aleatoric uncertainty. ' ...
            'Increase the mc_fit number_of_tries option by a factor of ~%d, or otherwise increase the central-model sample density.\n\n' ...
            'Continue anyway?'], ...
            n_pre, n_req, 100*opts.target_rse, fpre);
        do_continue = mc_fit_dialogs('ask_small_sample_continue', msg);
        if ~do_continue
            error('mc_fit:InsufficientCentralSampleBeforeFiltering', ...
                'Aborted because too few central models were available before imprecision filtering.');
        end
    end

    retained_fraction_threshold = 0.05;
    retained_fraction = n_post / max(n_pre,1);
    retained_fraction_is_small = retained_fraction < retained_fraction_threshold;

    if n_post < 3
        stats.aleatoric_equated_to_data = true;
        stats.aleatoric_covariance_available = false;

        if retained_fraction_is_small
            msg = sprintf(['Only %d central model(s) were retained after imprecision filtering; %d of %d available central model(s) were removed (%.1f%% retained).\n\n' ...
                'The retained central-model sample is too small to calculate an aleatoric covariance ellipse. ' ...
                'Because only %.1f%% of the available central models were retained, the filtering result suggests that the aleatoric uncertainty may be negligible relative to the data uncertainty. ' ...
                'This inference is based on the small proportion of retained models rather than on the absolute number of retained models alone.\n\n' ...
                'If desired, quantify the aleatoric uncertainty by re-running mc_fit with more central models (increase the number_of_tries option), ' ...
                'or by restricting the range of initial conditions around the current best central model.\n\n' ...
                'This case should be distinguished from a run in which few central models were available prior to imprecision filtering; here %d central model(s) were available before filtering.\n\n' ...
                'Continuing from the current condition will formally equate the total uncertainty ellipse to the data uncertainty ellipse and draw it that way.'], ...
                n_post, n_removed, n_pre, 100*retained_fraction, 100*retained_fraction, n_pre);
        else
            msg = sprintf(['Only %d central model(s) were retained after imprecision filtering; %d of %d available central model(s) were removed (%.1f%% retained).\n\n' ...
                'The retained central-model sample is too small to calculate an aleatoric covariance ellipse. Because %.1f%% of the available central models were retained, ' ...
                'this condition should not be interpreted as evidence that the aleatoric uncertainty is negligible relative to the data uncertainty; it may instead reflect inadequate sampling of the central-model space near the imprecision-defined acceptance interval.\n\n' ...
                'If desired, quantify the aleatoric uncertainty by re-running mc_fit with more central models (increase the number_of_tries option), ' ...
                'or by restricting the range of initial conditions around the current best central model.\n\n' ...
                'Continuing from the current condition will formally equate the total uncertainty ellipse to the data uncertainty ellipse and draw it that way, because no aleatoric covariance ellipse can be estimated from fewer than three retained central models.'], ...
                n_post, n_removed, n_pre, 100*retained_fraction, 100*retained_fraction);
        end

        warning('mc_fit:AleatoricCovarianceNotEstimable', '%s', msg);
        try
            choice = questdlg(msg, 'MC_fit aleatoric covariance warning', ...
                'Continue','Abort','Continue');
            if isempty(choice) || strcmp(choice,'Abort')
                error('mc_fit:UserAbort','User aborted after aleatoric covariance warning.');
            end
        catch ME
            if strcmp(ME.identifier,'mc_fit:UserAbort')
                rethrow(ME);
            end
            resp = strtrim(input(sprintf('\n%s\n\nContinue? [Y]/N: ',msg),'s'));
            if strcmpi(resp,'N') || strcmpi(resp,'NO')
                error('mc_fit:UserAbort','User aborted after aleatoric covariance warning.');
            end
        end
    elseif n_post < n_req
        fpost = ceil(n_req / max(n_post,1));

        if retained_fraction_is_small
            msg = sprintf(['Only %d central model(s) were retained after imprecision filtering; %d of %d available central model(s) were removed (%.1f%% retained). ' ...
                'At least %d retained models are required for the estimated aleatoric standard deviation to be accurate to within %.0f%%.\n\n' ...
                'The aleatoric covariance ellipse can still be calculated, but it is poorly constrained. Because only %.1f%% of the available central models were retained, the filtering result suggests that the aleatoric uncertainty may be negligible relative to the data uncertainty, even though the covariance estimate itself is imprecise.\n\n' ...
                'If desired, improve the aleatoric uncertainty estimate by re-running mc_fit with more central models (increase number_of_tries by a factor of ~%d), ' ...
                'or by restricting the range of initial conditions around the current best central model.\n\n' ...
                'Continue anyway?'], ...
                n_post, n_removed, n_pre, 100*retained_fraction, n_req, 100*opts.target_rse, 100*retained_fraction, fpost);
        else
            msg = sprintf(['Only %d central model(s) were retained after imprecision filtering; %d of %d available central model(s) were removed (%.1f%% retained). ' ...
                'At least %d retained models are required for the estimated aleatoric standard deviation to be accurate to within %.0f%%.\n\n' ...
                'The aleatoric covariance ellipse can still be calculated, but it is poorly constrained. Because %.1f%% of the available central models were retained, ' ...
                'this condition should not be interpreted as evidence that the aleatoric uncertainty is negligible relative to the data uncertainty; it may instead reflect inadequate sampling of the central-model space near the imprecision-defined acceptance interval.\n\n' ...
                'Increase number_of_tries by a factor of ~%d, or restrict the range of initial conditions around the current best central model to improve sampling near the optimum.\n\n' ...
                'Continue anyway?'], ...
                n_post, n_removed, n_pre, 100*retained_fraction, n_req, 100*opts.target_rse, 100*retained_fraction, fpost);
        end

        do_continue = mc_fit_dialogs('ask_small_sample_continue', msg);

        if ~do_continue
            error('mc_fit:InsufficientCentralSampleAfterFiltering', ...
                'Aborted because too few central models were retained after imprecision filtering.');
        end

    end
end

if ~handled_sample_warning && n_cov_points > 0
    f = ceil(n_req / max(n_cov_points,1));

    if n_cov_points < n_req
        if model == "Perturbation"
            knob = 'number_of_perturbations';
        else
            knob = 'number_of_tries';
        end

        msg = sprintf(['Only %d points are available for the variance estimates; ' ...
            'at least %d are required for the estimated standard deviation to be accurate to within %.0f%%.\n\n' ...
            'Increase %s by a factor of ~%d to obtain sufficient sampling.\n\n' ...
            'Continue anyway?'], ...
            n_cov_points, n_req, 100*opts.target_rse, knob, f);

        do_continue = mc_fit_dialogs('ask_small_sample_continue', msg);
        if ~do_continue
            error('mc_fit:InsufficientCovarianceSample', ...
                'Aborted because the covariance sample is too small.');
        end
    end
end

if n_cov_points >= 3
    data_xy = [x(idx_sigma_cov), y(idx_sigma_cov)];
    if filter && model == "Central" && ~isempty(pcovar)
        [handles.cov, handles.cov_total, ~, pcovar, handles.data_cov] = ...
            function_to_plot_covariance_ellipse(data_xy, stats.best_z, 3, dlw, pcovar, true, opts.sigma_level, [stats.best_x stats.best_y]);
    else
        jplot = 3;
        if model == "Perturbation", jplot = 2; end
        [handles.cov, handles.cov_total, ~, pcovar, handles.data_cov] = ...
            function_to_plot_covariance_ellipse(data_xy, stats.best_z, jplot, dlw, pcovar, false, opts.sigma_level, [stats.best_x stats.best_y]);
    end
elseif filter && model == "Central" && ~isempty(pcovar)
    % The aleatoric covariance ellipse cannot be calculated from fewer than
    % three retained central models.  In the negligible-aleatoric limiting case,
    % draw the total uncertainty ellipse as the perturbation-derived data
    % uncertainty ellipse.
    stats.aleatoric_equated_to_data = true;
    stats.aleatoric_covariance_available = false;
    handles.data_cov = plot_covariance_matrix_ellipse(h(1), pcovar, stats.best_x, stats.best_y, stats.best_z, '-r', dlw, opts.sigma_level);
    handles.cov_total = plot_covariance_matrix_ellipse(h(1), pcovar, stats.best_x, stats.best_y, stats.best_z, '-k', dlw, opts.sigma_level);
end

if ~isfield(opts,'plotcov') || ~opts.plotcov
    local_hide_if_graphics(handles.cov)
    local_hide_if_graphics(handles.cov_total)
    local_hide_if_graphics(handles.data_cov)
end

% MPP
if opts.plotMPP && model == "Central" && ~filter && all(~isnan([MPPX MPPY MPPZ]))
    handles.mpp = plot_mpp(h(1), MPPX, MPPY, MPPZ);
end

% Legend
legend_entries = build_legend_entries(model, filter, stats, sname, uname, zname, ...
    h_dum, handles, opts.printz, opts.print, opts);
hlist = legend_entries.hlist;
tlist = legend_entries.tlist;
if ~isempty(hlist)
    hleg = legend(h(1), hlist, tlist);
    hleg.Location = 'northwest';
end

title(h(1), make_title_text(model,filter,opts))

h(1).XAxis.Exponent = 0;
h(1).YAxis.Exponent = 0;
h(1).ZAxis.Exponent = 0;

end


function local_hide_if_graphics(h)
    if isempty(h), return, end
    if all(isgraphics(h))
        set(h, 'Visible', 'off')
    end
end

function hcov = plot_covariance_matrix_ellipse(ax, C, X0, Y0, Z0, line_spec, lw, sigma_level)
    if isempty(C) || any(~isfinite(C(:)))
        hcov = gobjects(1,1);
        return
    end

    C = (C + C.')/2;
    [V,D] = eig(C);
    d = max(diag(D), 0);
    [largest_eigenval, idx] = max(d);
    largest_eigenvec = V(:,idx);
    smallest_eigenval = min(d);

    phi = atan2(largest_eigenvec(2), largest_eigenvec(1));
    if phi < 0
        phi = phi + 2*pi;
    end

    p = erf(sigma_level / sqrt(2));
    c = -2 * log(1 - p);
    chisquare_val = sqrt(c);
    theta = linspace(0,2*pi);

    a = chisquare_val * sqrt(largest_eigenval);
    b = chisquare_val * sqrt(smallest_eigenval);
    xr = a * cos(theta);
    yr = b * sin(theta);
    R = [cos(phi) sin(phi); -sin(phi) cos(phi)];
    rell = [xr; yr]' * R;
    zel = Z0 * ones(size(rell,1),1);

    hcov = plot3(ax, rell(:,1)+X0, rell(:,2)+Y0, zel, line_spec, 'LineWidth', lw);
end


function h0 = plot_best_model(ax,x,y,z,fit,idx,dlw,filled_symbols_for_nofit)
    h0 = plot_best_model_coords(ax, x(idx), y(idx), z(idx), fit(idx), dlw, filled_symbols_for_nofit);
end

function h0 = plot_best_model_coords(ax,x0,y0,z0,fit0,dlw,filled_symbols_for_nofit)
    h0 = plot_model_marker(ax, x0, y0, z0, fit0, "#EDB120", 'o', 25, dlw, filled_symbols_for_nofit);
end

function [hdata, stats] = draw_perturbation_errorbars(ax, x, y, z, idx, stats, dlw, sigma_level, cap_fraction)
    stats.sigx = std(x(idx), 0, 'omitnan');
    stats.sigy = std(y(idx), 0, 'omitnan');
    stats.epsz = std(z(idx), 0, 'omitnan');
    stats.filter_half_width_z = stats.epsz;
    stats.pcovar = [];
    % Center perturbation error bars on the best central model when supplied
    % by mc_fit_plot's redraw pass; otherwise fall back to the perturbation mean.
    if isfield(stats,'best_x') && isfinite(stats.best_x) && ...
            isfield(stats,'best_y') && isfinite(stats.best_y) && ...
            isfield(stats,'best_z') && isfinite(stats.best_z)
        ctrx = stats.best_x;
        ctry = stats.best_y;
        ctrz = stats.best_z;
    else
        ctrx = mean(x(idx), 'omitnan');
        ctry = mean(y(idx), 'omitnan');
        ctrz = mean(z(idx), 'omitnan');
    end
    stats.sigx = sigma_level * stats.sigx;
    stats.sigy = sigma_level * stats.sigy;
    stats.epsz = sigma_level * stats.epsz;
    hdata = function_errorbar3_caps(ax, ctrx, ctry, ctrz, ...
        stats.sigx, stats.sigy, stats.epsz, 'Color', 'red', 'LineWidth', dlw, ...
        'cap_fraction', cap_fraction);
end

function [himp, htot, stats] = draw_central_errorbars(ax, x, y, z, idx, idx_best, pcovar, filter, stats, dlw, sigma_level, cap_fraction, central_epsz_override)
    data_xyz = [x(idx), y(idx), z(idx)];
    n_finite = sum(all(isfinite(data_xyz), 2));

    stats.sigx = std(x(idx), 0, 'omitnan');
    stats.sigy = std(y(idx), 0, 'omitnan');
    stats.epsz = std(z(idx), 0, 'omitnan');
    use_epsz_override = nargin >= 13 && isfinite(central_epsz_override);
    if use_epsz_override
        % Some filtered central analyses define the accepted z-range by the
        % perturbation-derived misfit imprecision, but the retained central
        % sample may not provide a meaningful independent aleatoric z scatter.
        % In that case, use the supplied perturbation imprecision as the
        % aleatoric z component so the vertical error bars and legend remain
        % consistent with the plotted acceptance interval.
        stats.epsz = central_epsz_override / sigma_level;
    end

    if isempty(pcovar)
        pert_sigx = 0;
        pert_sigy = 0;
        stats.pcovar = [];
    else
        pert_sigx = sqrt(max(pcovar(1,1),0));
        pert_sigy = sqrt(max(pcovar(2,2),0));
        stats.pcovar = pcovar;
    end

    no_aleatoric_covariance = filter && n_finite < 3 && ~isempty(pcovar);
    if no_aleatoric_covariance
        % Aleatoric covariance ellipse cannot be estimated from fewer than three
        % retained central models.  Use the documented negligible-aleatoric
        % limiting case: total uncertainty equals the perturbation-derived
        % data uncertainty.  Keep a blue zero-valued aleatoric legend entry.
        stats.sigx = 0;
        stats.sigy = 0;
        if ~use_epsz_override
            stats.epsz = 0;
        end
        stats.aleatoric_equated_to_data = true;
        stats.aleatoric_covariance_available = false;
    end

    stats.sigxt = sqrt(stats.sigx.^2 + pert_sigx.^2);
    stats.sigyt = sqrt(stats.sigy.^2 + pert_sigy.^2);

    stats.sigx = sigma_level * stats.sigx;
    stats.sigy = sigma_level * stats.sigy;
    stats.epsz = sigma_level * stats.epsz;
    stats.sigxt = sigma_level * stats.sigxt;
    stats.sigyt = sigma_level * stats.sigyt;

    if filter
        real_epsz = max(0, stats.epsz);
        htot = function_errorbar3_caps(ax, x(idx_best), y(idx_best), z(idx_best), ...
            stats.sigxt, stats.sigyt, real_epsz, 'Color', 'k', 'LineWidth', dlw, ...
            'cap_fraction', cap_fraction);
    else
        htot = [];
    end

    if no_aleatoric_covariance
        % No aleatoric covariance/error bar can be estimated.  Create an
        % explicit blue legend proxy (not a plotted ellipse/error bar) so
        % the filtered central-analysis legend still reports
        % sigma_XY,aleatoric = 0.  A proxy is used instead of
        % function_errorbar3_caps because those internal line objects have
        % HandleVisibility off and may be omitted by legend.
        % Use a zero-length visible blue line at the best-model coordinate
        % as the legend handle.  NaN-only proxy lines can be omitted from
        % legends in some MATLAB graphics states, whereas a zero-length line
        % gives the legend a real handle without visibly changing the plot.
        xb = x(idx_best); yb = y(idx_best); zb = z(idx_best);
        himp = plot3(ax, [xb xb], [yb yb], [zb zb], '-b', 'LineWidth', dlw, ...
            'DisplayName', 'Aleatoric uncertainty');
    else
        himp = function_errorbar3_caps(ax, x(idx_best), y(idx_best), z(idx_best), ...
            stats.sigx, stats.sigy, real_epsz, 'Color', 'blue', 'LineWidth', dlw, ...
            'cap_fraction', cap_fraction);
    end
end

function stats = compute_stats_only(x, y, z, idx, stats, model, filter, idx_best, pcovar)
    stats.sigx = std(x(idx), 0, 'omitnan');
    stats.sigy = std(y(idx), 0, 'omitnan');
    stats.epsz = std(z(idx), 0, 'omitnan');
    if model == "Perturbation"
        stats.filter_half_width_z = stats.epsz;
        stats.pcovar = [];
    elseif filter
        if isempty(pcovar)
            pert_sigx = 0;
            pert_sigy = 0;
            stats.pcovar = [];
        else
            pert_sigx = sqrt(max(pcovar(1,1),0));
            pert_sigy = sqrt(max(pcovar(2,2),0));
            stats.pcovar = pcovar;
        end
        stats.sigxt = sqrt(stats.sigx.^2 + pert_sigx.^2);
        stats.sigyt = sqrt(stats.sigy.^2 + pert_sigy.^2);
    end
    stats.best_x = x(idx_best);
    stats.best_y = y(idx_best);
    stats.best_z = z(idx_best);
end

function [h1a,h1b,h_dum] = plot_trials(ax,x,y,z,fit,ind,ec,mt,dlw,filled_symbols_for_nofit)
    h1a = gobjects(1,1);
    h1b = gobjects(1,1);
    h_dum = plot(ax,nan,nan,'LineStyle','none');

    sz = 6;
    lw = 0.5*dlw;

    jnd = ind(fit(ind)==1);
    if ~isempty(jnd)
        h1a = plot_model_marker(ax, x(jnd), y(jnd), z(jnd), 1, ec, mt, sz, lw, filled_symbols_for_nofit, ec);
        h_dum = h1a;
    end

    jnd = ind(fit(ind)==0);
    if ~isempty(jnd)
        h1b = plot_model_marker(ax, x(jnd), y(jnd), z(jnd), 0, ec, mt, sz, lw, filled_symbols_for_nofit, ec);
        if ~isgraphics(h_dum)
            h_dum = h1b;
        end
    end
end

function h2 = plot_best_overall(ax,x,y,z,fit,idx,dlw,filled_symbols_for_nofit)
    h2 = plot_model_marker(ax, x(idx), y(idx), z(idx), fit(idx), "#0072BD", 'o', 25, dlw, filled_symbols_for_nofit);
end

function h = plot_model_marker(ax,x,y,z,fit_value,fill_color,marker,marker_size,line_width,filled_symbols_for_nofit,edge_color)
    if nargin < 11
        edge_color = 'k';
    end
    [fc, ec] = marker_fill_from_fit(fit_value, fill_color, edge_color, filled_symbols_for_nofit);
    h = plot_marker(ax, x, y, z, marker, fc, ec, marker_size, line_width);
end

function h = plot_marker(ax,x,y,z,marker,face_color,edge_color,marker_size,line_width)
    h = plot3(ax, x, y, z, marker, ...
        'MarkerFaceColor', face_color, ...
        'MarkerEdgeColor', edge_color, ...
        'MarkerSize', marker_size, ...
        'LineWidth', line_width);
end

function [fc, ec] = marker_fill_from_fit(fit_value, fill_color, edge_color, filled_symbols_for_nofit)
    % Map the filled/open marker state to the data-fit flag.  By default,
    % filled symbols identify models that do not fit all analytical data.
    is_fit = (fit_value == 1);
    make_filled = xor(is_fit, filled_symbols_for_nofit);
    if make_filled
        fc = fill_color;
        ec = edge_color;
    else
        fc = 'none';
        ec = edge_color;
    end
end

function hmpp = plot_mpp(ax,MPPX,MPPY,MPPZ)
    hmpp = plot_marker(ax, MPPX, MPPY, MPPZ, 'd', 'k', 'k', 15, 0.5);
end

function ttl = make_title_text(model,filter,opts)
    if model == "Central" && filter
        ttl = ['Filtered ' char(model) ' Analysis'];
    else
        ttl = [char(model) ' Analysis'];
    end

    if isstruct(opts) && isfield(opts,'title_prefix') && strlength(string(opts.title_prefix)) > 0
        ttl = char(string(opts.title_prefix) + ": " + string(ttl));
    end
end

function format_axes(ax,xname,yname,zname,logz)
    view(ax,0,90)
    xlabel(ax,xname)
    ylabel(ax,yname)
    zlabel(ax,zname)
    if logz
        zscale(ax,'log')
    end
    axis(ax,'square')
    box(ax,'on')
end


function entries = build_legend_entries(model, filter, stats, sname, uname, zname, h_dum, handles, printz, printflag, opts)
    entries.hlist = first_graphics_handle(handles.best);
    entries.tlist = {make_best_model_text(stats, sname, uname, zname, printz, opts)};

    if model == "Central" && filter
        if isfield(opts,'plotcov') && opts.plotcov
            h_data = first_graphics_handle(handles.data_cov);
        else
            h_data = gobjects(1,1);
        end
        if ~isgraphics(h_data)
            h_data = first_graphics_handle(handles.data);
        end
        if isgraphics(h_data)
            entries.hlist = [entries.hlist h_data];
            entries.tlist = [entries.tlist {make_component_text("data", stats, sname, uname, zname, printz, printflag, opts)}];
        end
        %eliminate best central and data uncertainty from filtered result:
        entries.hlist = [];
        entries.tlist = [];
        % Imprecision is keyed by the blue covariance ellipse/error bar when
        % present.  In the negligible-aleatoric limiting case the blue handle
        % is a dummy object so the legend still reports sigma_XY,aleatoric = 0.
        if isfield(opts,'plotcov') && opts.plotcov
            h_imp = first_graphics_handle(handles.cov);
        else
            h_imp = gobjects(1,1);
        end
        if ~isgraphics(h_imp)
            h_imp = first_graphics_handle(handles.imprecision);
        end
        if ~isgraphics(h_imp) && isfield(stats,'aleatoric_equated_to_data') && stats.aleatoric_equated_to_data
            % Last-resort legend proxy for the negligible-aleatoric limiting
            % case.  This prevents the aleatoric legend entry from being lost
            % if no blue covariance/error-bar handle was created upstream.
            h_best_for_ax = first_graphics_handle(handles.best);
            if isgraphics(h_best_for_ax)
                ax_proxy = ancestor(h_best_for_ax,'axes');
            else
                ax_proxy = gca;
            end
            h_imp = plot3(ax_proxy, [stats.best_x stats.best_x], [stats.best_y stats.best_y], ...
                [stats.best_z stats.best_z], '-b', 'LineWidth', get(groot,'DefaultLineLineWidth'), ...
                'DisplayName', 'Aleatoric uncertainty');
        end
        if isgraphics(h_imp)
            entries.hlist = [entries.hlist h_imp];
            entries.tlist = [entries.tlist {make_component_text("imprecision", stats, sname, uname, zname, printz, printflag, opts)}];
        end

        if isfield(opts,'plotcov') && opts.plotcov
            h_tot = first_graphics_handle(handles.cov_total);
        else
            h_tot = gobjects(1,1);
        end
        if ~isgraphics(h_tot)
            h_tot = first_graphics_handle(handles.total);
        end
        if isgraphics(h_tot)
            entries.hlist = [entries.hlist h_tot];
            entries.tlist = [entries.tlist {make_component_text("total", stats, sname, uname, zname, printz, printflag, opts)}];
        end
    else
        trial_handle = first_graphics_handle(h_dum);
        if model == "Perturbation"
            if isfield(opts,'plotcov') && opts.plotcov
                h_cov = first_graphics_handle(handles.cov);
                if isgraphics(h_cov)
                    trial_handle = h_cov;
                end
            end
        end
        if isgraphics(trial_handle)
            entries.hlist = [entries.hlist trial_handle];
            entries.tlist = [entries.tlist {make_trial_text(model, filter, stats, sname, uname, zname, printz, printflag, opts)}];

            if model == "Perturbation" && isfield(stats,'epsz') && isfinite(stats.epsz)
                % Add the imprecision value as a text-only legend entry. Do
                % not reuse the trial marker handle here: when some models
                % fit all data, h_dum may point to a real model marker, which
                % incorrectly puts that marker next to the imprecision text.
                h_imp_text = make_blank_legend_handle(ancestor(trial_handle,'axes'));
                if isgraphics(h_imp_text)
                    entries.hlist = [entries.hlist h_imp_text];
                    entries.tlist = [entries.tlist {make_imprecision_z_text(stats, zname, opts)}];
                else
                    % Fallback: append imprecision to the data-uncertainty entry.
                    entries.tlist{end} = [entries.tlist{end} newline make_imprecision_z_text(stats, zname, opts)];
                end
            end
        end
    end

    h_overall = first_graphics_handle(handles.overall);
    if isgraphics(h_overall)
        entries.hlist = [entries.hlist h_overall];
        entries.tlist = [entries.tlist {make_best_overall_text(stats, sname, uname, zname, printz, opts)}];
    end

    h_mpp = first_graphics_handle(handles.mpp);
    if isgraphics(h_mpp)
        entries.hlist = [entries.hlist h_mpp];
        entries.tlist = [entries.tlist {make_mpp_text(sname, uname, h_mpp.XData, h_mpp.YData, opts)}];
    end
end

function h = first_graphics_handle(hin)
    h = gobjects(1,1);
    if isempty(hin)
        return
    end

    if isstruct(hin)
        fn = fieldnames(hin);
        for jf = 1:numel(fn)
            h = first_graphics_handle(hin.(fn{jf}));
            if isgraphics(h)
                return
            end
        end
        return
    end

    try
        hin = hin(:);
        k = find(isgraphics(hin), 1, 'first');
        if ~isempty(k)
            h = hin(k);
        end
    catch
        h = gobjects(1,1);
    end
end



function h = make_blank_legend_handle(ax)
    if nargin < 1 || ~isgraphics(ax,'axes')
        ax = gca;
    end
    h = plot3(ax, nan, nan, nan, ...
        'LineStyle','none', 'Marker','none', 'Color','none', ...
        'HandleVisibility','on');
end

function strg = make_imprecision_z_text(stats, zname, opts)
    %#ok<INUSD>  % zname retained for interface consistency with older calls.
    strg = ['\bfImprecision\rm' newline ...
        '\rm\epsilon_{ln(Misfit)} = ' sigfig(stats.epsz,opts.legend_sigma_digits)];
end

function strg = make_best_model_text(stats, sname, uname, zname, printz, opts)
    xlab = legend_var_label(sname(1));
    ylab = legend_var_label(sname(2));
    xunit = legend_plain_label(uname(1));
    yunit = legend_plain_label(uname(2));
    zlab = legend_plain_label(zname);

    strg = ['\bfBest central model\rm' newline ...
        '\rm' xlab ' = ' value_with_unit(stats.best_x, xunit, opts.legend_coord_digits) newline ...
        '\rm' ylab ' = ' value_with_unit(stats.best_y, yunit, opts.legend_coord_digits)];
    if printz
        strg = [strg newline '\rm' zlab ' = ' sigfig(stats.best_z,opts.legend_coord_digits)];
    end
end

function strg = make_trial_text(model, filter, stats, sname, uname, zname, printz, printflag, opts)
    if model == "Central" && ~filter
        strg = make_component_text("tries", stats, sname, uname, zname, printz, printflag, opts);
    else
        strg = make_component_text("perturbation", stats, sname, uname, zname, printz, printflag, opts);
    end
end

function strg = make_component_text(kind, stats, sname, uname, zname, printz, printflag, opts)
    xlab = legend_var_label(sname(1));
    ylab = legend_var_label(sname(2));
    xunit = legend_plain_label(uname(1));
    yunit = legend_plain_label(uname(2));
    zlab = legend_plain_label(zname);

    switch string(kind)
        case "tries"
            head = '\bfTries\rm';
            xval = stats.sigx; yval = stats.sigy; zval = stats.epsz;
            xtag = 'tries'; ytag = 'tries'; ztag = 'tries'; zsym = '\epsilon';
        case "perturbation"
            if printflag
                head = '\bfData uncertainty\rm';
            else
                head = '\bfData\rm';
            end
            % For perturbation plots, the x-y scatter defines the data-
            % uncertainty ellipse.  The misfit scatter (epsz) is reported
            % separately as imprecision, so do not repeat it here.
            xval = stats.sigx; yval = stats.sigy; zval = [];
            xtag = 'data'; ytag = 'data'; ztag = 'data'; zsym = '\epsilon';
        case "data"
            if printflag
                head = '\bfData uncertainty\rm';
            else
                head = '\bfData\rm';
            end
            xval = sqrt(max(pcovar_diag(stats,'x'),0));
            yval = sqrt(max(pcovar_diag(stats,'y'),0));
            zval = []; xtag = 'data'; ytag = 'data'; ztag = 'data'; zsym = '\epsilon';
        case "imprecision"
            head = '\bfAleatoric uncertainty\rm';
            xval = stats.sigx; yval = stats.sigy; zval = stats.epsz;
            xtag = 'aleatoric'; ytag = 'aleatoric'; ztag = 'aleatoric'; zsym = '\epsilon';
        case "total"
            if printflag
                head = '\bfTotal uncertainty\rm';
            else
                head = '\bfTotal\rm';
            end
            xval = stats.sigxt; yval = stats.sigyt; zval = [];
            xtag = 'total'; ytag = 'total'; ztag = 'total'; zsym = '\epsilon';
        otherwise
            head = '\bfComponent\rm';
            xval = stats.sigx; yval = stats.sigy; zval = [];
            xtag = 'component'; ytag = 'component'; ztag = 'component'; zsym = '\epsilon';
    end

    strg = [head newline ...
        '\rm\sigma_{' component_subscript(xlab, xtag) '} = ' value_with_unit(xval, xunit, opts.legend_sigma_digits) newline ...
        '\rm\sigma_{' component_subscript(ylab, ytag) '} = ' value_with_unit(yval, yunit, opts.legend_sigma_digits)];

    if printz && ~isempty(zval)
        strg = [strg newline '\rm' zsym '_{' component_subscript(zlab, ztag) '} = ' sigfig(zval,opts.legend_sigma_digits)];
    end
end

function val = pcovar_diag(stats, which)
    if isfield(stats, 'pcovar') && ~isempty(stats.pcovar)
        if which == 'x'
            val = stats.pcovar(1,1);
        else
            val = stats.pcovar(2,2);
        end
    else
        val = 0;
    end
end

function strg = make_best_overall_text(stats, sname, uname, zname, printz, opts)
    xlab = legend_var_label(sname(1));
    ylab = legend_var_label(sname(2));
    xunit = legend_plain_label(uname(1));
    yunit = legend_plain_label(uname(2));
    zlab = legend_plain_label(zname);

    strg = ['\bfBest overall model\rm' newline ...
        '\rm' xlab ' = ' value_with_unit(stats.overall_x, xunit, opts.legend_coord_digits) newline ...
        '\rm' ylab ' = ' value_with_unit(stats.overall_y, yunit, opts.legend_coord_digits)];
    if printz
        strg = [strg newline '\rm' zlab ' = ' sigfig(stats.overall_z,opts.legend_coord_digits)];
    end
end
function strg = make_mpp_text(sname, uname, MPPX, MPPY, opts)
    xlab = legend_var_label(sname(1));
    ylab = legend_var_label(sname(2));
    xunit = legend_plain_label(uname(1));
    yunit = legend_plain_label(uname(2));

    strg = ['\bfMPP\rm' newline ...
        '\rm' xlab ' = ' value_with_unit(MPPX, xunit, opts.legend_coord_digits) newline ...
        '\rm' ylab ' = ' value_with_unit(MPPY, yunit, opts.legend_coord_digits)];
end

function out = component_subscript(label, tag)
    out = [legend_subscript_label(label) ',' char(string(tag))];
end

function out = legend_subscript_label(in)
    base = legend_plain_label(in);
    if startsWith(base, 'ln(') || startsWith(base, 'log(')
        out = base;
    elseif contains(base, '_')
        parts = split(string(base), '_');
        head = char(parts(1));
        tail = char(join(parts(2:end), '_'));
        out = ['\it' head '_{\rm' tail '}'];
    else
        out = ['\it' base '\rm'];
    end
end

function s = sigfig(val, n)
    if nargin < 2 || isempty(n)
        n = 4;
    end
    if isempty(val) || ~isfinite(val)
        s = num2str(val);
    else
        s = sprintf('%.*g', n, val);
    end
end

function out = value_with_unit(val, unit, nfmt)
    if nargin < 3 || isempty(nfmt)
        nfmt = 4;
    end
    sval = sigfig(val, nfmt);
    u = legend_plain_label(unit);
    if isempty(u)
        out = sval;
    else
        out = [sval ' ' u];
    end
end

function out = legend_plain_label(in)
    out = char(string(in));
    out = regexprep(out, '\\(it|rm|bf)\s*', '');
    out = strtrim(out);
end

function out = legend_var_label(in)
    base = legend_plain_label(in);
    if contains(base, '_')
        parts = split(string(base), '_');
        head = char(parts(1));
        tail = char(join(parts(2:end), '_'));
        out = ['\it' head '_{\rm' tail '}'];
    else
        out = ['\it' base '\rm'];
    end
end
