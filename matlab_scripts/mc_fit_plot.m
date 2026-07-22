% Script to plot mc_fit *.pts format files. The script prompts
% for the *_perturbed.pts file and automatically reads the corresponding
% *_central.pts file.
%
% JADC March 20, 2026
% Refactored by ChatGPT April 5, 2026

% -------------------------------------------------------------------------
% switches
print   = true;    % false: projection/PPT settings; true: print settings
outPDF  = false;   % false: png bitmap; true: vector pdf
plotMPP = false;   % plot MPP in central model
plotcov = true;    % plot covariance ellipse(s)
ploterr = true;    % plot error bars
plotbst = true;    % plot best overall model if better than best central
auto_filled_symbols_for_nofit = true;  % true: infer marker convention from the fit-all-data distribution
filled_symbols_for_nofit = false;      % used only when auto_filled_symbols_for_nofit is false; true fills non-fitting models
printz  = false;   % print misfit/bayes scores in legends
maksub  = true;    % make subplots
fitonly_prompt = false;   % prompt whether to retain only fitting models
manual_filter  = false;   % use manual filter dialog
filter         = true;    % automatically filter central models to within perturbation-derived imprecision
% -------------------------------------------------------------------------
%                           unit conversions
convert_to_MPa  = false;
convert_to_kbar = true;
convert_to_GPa  = false;
convert_to_C    = false;
% -------------------------------------------------------------------------
% user parameters
sigma_level    = 1;      % 1 -> 68%-style scaling, 2 -> 95%-style scaling
target_rse    = 0.20;   % relative precision threshold for std estimate warnings; central samples are checked before and after imprecision filtering
cap_fraction   = 0.10;
legend_coord_digits = 4;   % significant digits for legend coordinates
legend_sigma_digits = 2;   % significant digits for legend uncertainties   % cap size relative to orthogonal error
outlier_filter = true;  % filter perturbation outliers before estimating imprecision
outlier_filter_iterations = 10; % maximum iterations for iterative outlier filtering
outlier_filter_convergence = 0.01; % stop if fractional change in retained points is below this value
outlier_filter_z = true;       % robust one-sided z filter on ln(Misfit)/ln(Bayes); rejects only the bad tail
outlier_filter_delta_z = false;  % reject models farther than outlier_filter_delta_z_max from the best z
outlier_filter_delta_z_max = 3.91; % maximum allowed e-fold distance from best ln(Misfit)/ln(Bayes)
% -------------------------------------------------------------------------
% hardening: catch missing or invalid user options early
assert(exist('sigma_level','var') == 1, 'sigma_level not defined')
assert(exist('target_rse','var') == 1, 'target_rse not defined')
assert(exist('cap_fraction','var') == 1, 'cap_fraction not defined')
assert(exist('outlier_filter','var') == 1, 'outlier_filter not defined')
assert(exist('outlier_filter_iterations','var') == 1, 'outlier_filter_iterations not defined')
assert(exist('outlier_filter_convergence','var') == 1, 'outlier_filter_convergence not defined')
assert(exist('outlier_filter_z','var') == 1, 'outlier_filter_z not defined')
assert(exist('outlier_filter_delta_z','var') == 1, 'outlier_filter_delta_z not defined')
assert(exist('outlier_filter_delta_z_max','var') == 1, 'outlier_filter_delta_z_max not defined')
assert(isscalar(sigma_level) && isnumeric(sigma_level) && isfinite(sigma_level) && sigma_level > 0, ...
    'sigma_level must be a finite positive scalar')
assert(isscalar(target_rse) && isnumeric(target_rse) && isfinite(target_rse) && target_rse > 0 && target_rse < 1, ...
    'target_rse must be a finite scalar between 0 and 1')
assert(isscalar(cap_fraction) && isnumeric(cap_fraction) && isfinite(cap_fraction) && cap_fraction >= 0, ...
    'cap_fraction must be a finite non-negative scalar')
assert(islogical(auto_filled_symbols_for_nofit) || (isscalar(auto_filled_symbols_for_nofit) && isnumeric(auto_filled_symbols_for_nofit)), ...
    'auto_filled_symbols_for_nofit must be a logical scalar')
assert(islogical(filled_symbols_for_nofit) || (isscalar(filled_symbols_for_nofit) && isnumeric(filled_symbols_for_nofit)), ...
    'filled_symbols_for_nofit must be a logical scalar')
assert(islogical(outlier_filter) || (isscalar(outlier_filter) && isnumeric(outlier_filter)), ...
    'outlier_filter must be a logical scalar')
assert(isscalar(outlier_filter_iterations) && isnumeric(outlier_filter_iterations) && ...
    isfinite(outlier_filter_iterations) && outlier_filter_iterations >= 0 && ...
    outlier_filter_iterations == floor(outlier_filter_iterations), ...
    'outlier_filter_iterations must be a non-negative integer')
assert(isscalar(outlier_filter_convergence) && isnumeric(outlier_filter_convergence) && ...
    isfinite(outlier_filter_convergence) && outlier_filter_convergence >= 0 && ...
    outlier_filter_convergence <= 1, ...
    'outlier_filter_convergence must be a finite scalar between 0 and 1')
assert(islogical(outlier_filter_z) || (isscalar(outlier_filter_z) && isnumeric(outlier_filter_z)), ...
    'outlier_filter_z must be a logical scalar')
assert(islogical(outlier_filter_delta_z) || (isscalar(outlier_filter_delta_z) && isnumeric(outlier_filter_delta_z)), ...
    'outlier_filter_delta_z must be a logical scalar')
assert(isscalar(outlier_filter_delta_z_max) && isnumeric(outlier_filter_delta_z_max) && ...
    isfinite(outlier_filter_delta_z_max) && outlier_filter_delta_z_max >= 0, ...
    'outlier_filter_delta_z_max must be a finite non-negative scalar')

% -------------------------------------------------------------------------
% initialization
first = true;
function_to_set_graphics_defaults(print)
figure(1); clf

xvar = [];
yvar = [];
zvar = [];
epsz = NaN;
xname = "";
yname = "";
zname = "";
sname = strings(1,3);
uname = strings(1,3);
file  = '';
path  = '';
pcovar = [];
limits = struct('minx',Inf,'maxx',-Inf,'miny',Inf,'maxy',-Inf,'minz',Inf,'maxz',-Inf);

opts = struct( ...
    'plotMPP', plotMPP, ...
    'plotcov', plotcov, ...
    'ploterr', ploterr, ...
    'plotbst', plotbst, ...
    'printz',  printz, ...
    'print',   print, ...
    'maksub',  maksub, ...
    'filled_symbols_for_nofit', logical(filled_symbols_for_nofit), ...
    'sigma_level', sigma_level, ...
    'legend_coord_digits', legend_coord_digits, ...
    'legend_sigma_digits', legend_sigma_digits, ...
    'cap_fraction', cap_fraction, ...
    'target_rse', target_rse, ...
    'outlier_filter', logical(outlier_filter), ...
    'outlier_filter_iterations', outlier_filter_iterations, ...
    'outlier_filter_convergence', outlier_filter_convergence, ...
    'outlier_filter_z', logical(outlier_filter_z), ...
    'outlier_filter_delta_z', logical(outlier_filter_delta_z), ...
    'outlier_filter_delta_z_max', outlier_filter_delta_z_max);

% -------------------------------------------------------------------------
% open the *_perturbed.pts file
[x,y,z,symb,fit,xname,yname,zname,sname,uname,~,xvar,yvar,zvar,file,path,model,filled_symbols_for_nofit_auto,~] = ...
    function_to_get_mc_fit_plot_file( ...
    '*_perturbed.pts', xvar, yvar, zvar, epsz, fitonly_prompt, ...
    manual_filter, filter, xname, yname, zname, sname, file, path, first, ...
    convert_to_MPa, convert_to_kbar, convert_to_GPa, convert_to_C, false, ...
    outlier_filter, outlier_filter_iterations, outlier_filter_convergence, ...
    outlier_filter_z, outlier_filter_delta_z, outlier_filter_delta_z_max);

file_root = regexprep(char(file), '(_central|_perturbed)\.pts$', '');
title_prefix = string(file_root);
opts.title_prefix = title_prefix;
if auto_filled_symbols_for_nofit
    opts.filled_symbols_for_nofit = logical(filled_symbols_for_nofit_auto);
end
opts_pert = opts;

% Optional fit-based filtering.
% In this workflow the perturbation question comes first. If accepted,
% the same data-fit filter is then applied to the central models without a
% second prompt. If declined, the central models are asked separately when
% their file is opened.
filter_by_fit = false;
fitonly_prompt_central = false;
pert_sym = 2 + double(zname == "ln(Bayes)") * 3;
idx_pert_trials = find(symb == pert_sym & ~isnan(x) & ~isnan(y) & ~isnan(z));
if ~isempty(idx_pert_trials)
    nfit = sum(fit(idx_pert_trials) == 1);
    ntrial = numel(idx_pert_trials);
    if nfit > 0 && nfit < ntrial
        filter_by_fit = mc_fit_dialogs('ask_fit_filter_single', nfit, ntrial, 'perturbation');
        if filter_by_fit
            idx = idx_pert_trials(fit(idx_pert_trials) ~= 1);
            x(idx) = NaN; y(idx) = NaN; z(idx) = NaN;
        else
            fitonly_prompt_central = true;
        end
    end
end

if filter
    subspec = {{1,2,1}};
else
    subspec = {{1,2,2}};
end

[plot1, pertStats, ~, pcovar, limits] = function_to_do_mc_fit_plots( ...
    x, y, z, symb, fit, xname, yname, zname, sname, uname, model, filter, ...
    pcovar, limits, opts, subspec);

% Cache perturbation data in case the central analysis revises the best
% central model after the "fit all data" choice. If that happens, redraw
% the perturbation subplot with the revised best-central marker/legend.
pert_cache = struct('x', x, 'y', y, 'z', z, 'symb', symb, 'fit', fit, ...
    'xname', xname, 'yname', yname, 'zname', zname, ...
    'sname', sname, 'uname', uname, 'model', model);
% -------------------------------------------------------------------------
% switch to central file
file = strrep(file, "perturbed", "central");
first = false;
fitonly_prompt = false;

[x,y,z,symb,fit,xname,yname,zname,sname,uname,~,xvar,yvar,zvar,file,path,model,filled_symbols_for_nofit_auto,imprecision_filter_info] = ...
    function_to_get_mc_fit_plot_file( ...
    '*_central.pts', xvar, yvar, zvar, pertStats.filter_half_width_z, ...
    fitonly_prompt_central, manual_filter, filter, xname, yname, zname, sname, file, path, first, ...
    convert_to_MPa, convert_to_kbar, convert_to_GPa, convert_to_C, filter_by_fit, ...
    outlier_filter, outlier_filter_iterations, outlier_filter_convergence, ...
    outlier_filter_z, outlier_filter_delta_z, outlier_filter_delta_z_max);

opts_central = opts;
opts_central.imprecision_filter_info = imprecision_filter_info;
if auto_filled_symbols_for_nofit
    opts_central.filled_symbols_for_nofit = logical(filled_symbols_for_nofit_auto);
end

if filter
    subspec = {{1,2,2}};
else
    subspec = {{1,2,1}};
end

[plot2, centralStats, centralHandles, pcovar, limits] = function_to_do_mc_fit_plots( ...
    x, y, z, symb, fit, xname, yname, zname, sname, uname, model, filter, ...
    pcovar, limits, opts_central, subspec);

% Redraw the perturbation subplot after the central analysis is finalized so
% that its best-central marker, legend, and covariance ellipse are centered on
% the final best central model. This matters if the fit-all-data choice changes
% the best central model, and is harmless when it does not.
opts_redraw = opts_pert;
opts_redraw.best_override = struct('x', centralStats.best_x, ...
    'y', centralStats.best_y, 'z', centralStats.best_z, ...
    'fit', centralStats.best_fit);
if ~isempty(plot1) && isgraphics(plot1(1))
    cla(plot1(1))
end
if filter
    redraw_subspec = {{1,2,1}};
else
    redraw_subspec = {{1,2,2}};
end
[plot1, pertStats, ~, pcovar, limits] = function_to_do_mc_fit_plots( ...
    pert_cache.x, pert_cache.y, pert_cache.z, pert_cache.symb, pert_cache.fit, ...
    pert_cache.xname, pert_cache.yname, pert_cache.zname, ...
    pert_cache.sname, pert_cache.uname, pert_cache.model, filter, ...
    pcovar, limits, opts_redraw, redraw_subspec);

% -------------------------------------------------------------------------
% overlay perturbation-derived error bars onto the central-model plot,
% centered on the best central model.
if ploterr && ~isempty(plot2) && isgraphics(plot2(1))
    hover = function_errorbar3_caps( ...
        plot2(1), centralStats.best_x, centralStats.best_y, centralStats.best_z, ...
        sigma_level * pertStats.sigx, sigma_level * pertStats.sigy, 0, ...
        'Color', 'red', 'LineWidth', get(groot,'DefaultLineLineWidth'), ...
        'cap_fraction', cap_fraction);

    % z is intentionally ignored for front/back ordering.
    function_stack_errorbar_groups_xy(centralHandles.imprecision, hover, ...
        centralStats.sigx, centralStats.sigy, pertStats.sigx, pertStats.sigy);
end

% -------------------------------------------------------------------------
% final touches
if ~all(isfinite([limits.minx limits.maxx limits.miny limits.maxy limits.minz limits.maxz]))
    valid = isfinite(x) & isfinite(y) & isfinite(z);
    if any(valid)
        limits.minx = min(x(valid)); limits.maxx = max(x(valid));
        limits.miny = min(y(valid)); limits.maxy = max(y(valid));
        limits.minz = min(z(valid)); limits.maxz = max(z(valid));
    end
end

apply_limits = mc_fit_dialogs('ask_apply_xyz_limits', 'Set the same X-Y-Z limits for all plots?');

if apply_limits
    xr0 = [limits.minx limits.maxx];
    yr0 = [limits.miny limits.maxy];
    zr0 = [limits.minz limits.maxz];
    [xr, yr, zr] = mc_fit_dialogs('ask_xyz_limits', xr0, yr0, zr0);

    if ~isempty(xr) && ~isempty(yr) && ~isempty(zr)
        xlim(plot1(1),xr); xlim(plot2(1),xr);
        ylim(plot1(1),yr); ylim(plot2(1),yr);
        zlim(plot1(1),zr); zlim(plot2(1),zr);
    end
end

file = strrep(file, "_central.pts", "");

mc_fit_export_plot(gcf, file, outPDF);

