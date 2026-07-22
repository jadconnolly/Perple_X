% Script to plot mc_fit *.pts files using the current shared plotting
% functions while preserving the original legacy workflow.
%
% Workflow:
%   1) Prompt for the *_central.pts file.
%   2) Optionally apply a manual z-range filter to the central results.
%   3) Automatically read the matching *_perturbed.pts file.
%   4) Plot central results first, then perturbation results.
%
% This replaces the former split implementation through
% function_for_mc_fit_plots_I.m and function_for_mc_fit_plots_II.m.
%
% JADC Aug 27, 2025
% Refactored April 11, 2026

clf

print = true;
outPDF = false;   % false: png bitmap; true: vector pdf
function_to_set_graphics_defaults(print)
figure(1); clf

% -------------------------------------------------------------------------
% switches
plotMPP = true;    % legacy script plotted MPP on the central XY panel
plotcov = true;
ploterr = false;   % legacy script did not draw error bars
plotbst = true;    % plot best overall perturbed model
auto_filled_symbols_for_nofit = true;  % true: infer marker convention from the fit-all-data distribution
filled_symbols_for_nofit = false;        % used only when auto_filled_symbols_for_nofit is false; true fills non-fitting models
printz  = true;
maksub  = true;

%                                       unit conversions
convert_to_MPa  = false;
convert_to_kbar = true;
convert_to_GPa  = false;
convert_to_C    = true;

sigma_level    = 1;
cap_fraction   = 0.10;
legend_coord_digits = 4;   % significant digits for legend coordinates
legend_sigma_digits = 2;   % significant digits for legend uncertainties
target_rse    = 0.20;
outlier_filter = true;    % used only when filter is true and model is perturbation
outlier_filter_iterations = 10;
outlier_filter_convergence = 0.01;
outlier_filter_z = true;       % robust one-sided z filter on ln(Misfit)/ln(Bayes); rejects only the bad tail
outlier_filter_delta_z = false;
outlier_filter_delta_z_max = 3.91;

% -------------------------------------------------------------------------
% hardening: catch missing or invalid user options early
assert(exist('sigma_level','var') == 1, 'sigma_level not defined')
assert(exist('target_rse','var') == 1, 'target_rse not defined')
assert(exist('cap_fraction','var') == 1, 'cap_fraction not defined')
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
assert(exist('outlier_filter','var') == 1, 'outlier_filter not defined')
assert(exist('outlier_filter_iterations','var') == 1, 'outlier_filter_iterations not defined')
assert(exist('outlier_filter_convergence','var') == 1, 'outlier_filter_convergence not defined')
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
assert(exist('outlier_filter_z','var') == 1, 'outlier_filter_z not defined')
assert(exist('outlier_filter_delta_z','var') == 1, 'outlier_filter_delta_z not defined')
assert(exist('outlier_filter_delta_z_max','var') == 1, 'outlier_filter_delta_z_max not defined')
assert(islogical(outlier_filter_z) || (isscalar(outlier_filter_z) && isnumeric(outlier_filter_z)), ...
    'outlier_filter_z must be a logical scalar')
assert(islogical(outlier_filter_delta_z) || (isscalar(outlier_filter_delta_z) && isnumeric(outlier_filter_delta_z)), ...
    'outlier_filter_delta_z must be a logical scalar')
assert(isscalar(outlier_filter_delta_z_max) && isnumeric(outlier_filter_delta_z_max) && ...
    isfinite(outlier_filter_delta_z_max) && outlier_filter_delta_z_max >= 0, ...
    'outlier_filter_delta_z_max must be a finite non-negative scalar')

% -------------------------------------------------------------------------
% initialization
xvar = 0; yvar = 0; zvar = 0;
file = ""; path = "";
xname = ""; yname = ""; zname = "";
sname = strings(1,3); uname = strings(1,3);
epsz = 0;
first = true;
fitonly_prompt = true;
manual_filter = true;   % preserve legacy manual filtering workflow
filter = false;         % no automatic perturbation-width filtering in legacy
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
% open the *_central.pts file first, with optional manual filtering
[x, y, z, symb, fit, xname, yname, zname, sname, uname, ~, xvar, yvar, zvar, ...
    file, path, model, filled_symbols_for_nofit_auto] = ...
    function_to_get_mc_fit_plot_file( ...
    '*_central.pts', xvar, yvar, zvar, epsz, fitonly_prompt, ...
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

subspec = {{2,5,[1,2,6,7]}};
[hxy, centralStats, centralHandles, pcovar, limits] = function_to_do_mc_fit_plots( ...
    x, y, z, symb, fit, xname, yname, zname, sname, uname, model, filter, ...
    pcovar, limits, opts, subspec);

% -------------------------------------------------------------------------
% create the legacy x-z and y-z side panels by copying the central XY plot
xz = subplot(2,5,3);
yz = subplot(2,5,8);

exclude_handles = [centralHandles.mpp; centralHandles.cov(:); centralHandles.cov_total(:)];
local_copy_plot_children(hxy(1), xz, exclude_handles)
local_copy_plot_children(hxy(1), yz, exclude_handles)

view(xz,0,0)
view(yz,90,0)

xlabel(xz,xname);
ylabel(xz,yname);
zlabel(xz,zname);
axis(xz,'square');
box(xz,'on');

xlabel(yz,xname);
ylabel(yz,yname);
zlabel(yz,zname);
axis(yz,'square');
box(yz,'on');

% -------------------------------------------------------------------------
% now open the matching perturbation file, reusing the selected variables
file = strrep(file, "central", "perturbed");
first = false;
fitonly_prompt = true;
manual_filter = false;

[x, y, z, symb, fit, ~, ~, ~, ~, ~, ~, ~, ~, ~, file, path, model, filled_symbols_for_nofit_auto] = ...
    function_to_get_mc_fit_plot_file( ...
    '*_perturbed.pts', xvar, yvar, zvar, epsz, fitonly_prompt, ...
    manual_filter, filter, xname, yname, zname, sname, file, path, first, ...
    convert_to_MPa, convert_to_kbar, convert_to_GPa, convert_to_C, false, ...
    outlier_filter, outlier_filter_iterations, outlier_filter_convergence, ...
    outlier_filter_z, outlier_filter_delta_z, outlier_filter_delta_z_max);

if auto_filled_symbols_for_nofit
    opts.filled_symbols_for_nofit = logical(filled_symbols_for_nofit_auto);
end

subspec = {{2,5,[4,5,9,10]}};
[hpert, pertStats, pertHandles, pcovar, limits] = function_to_do_mc_fit_plots( ...
    x, y, z, symb, fit, xname, yname, zname, sname, uname, model, filter, ...
    pcovar, limits, opts, subspec);

% -------------------------------------------------------------------------
% optional warning for too-small perturbation sample sizes

% Remove the automatically-added best-central legend entry from the
% perturbation panel, but keep the perturbation statistics and best-overall
% entry if present.
local_reset_perturbation_legend(hpert(1), pertHandles, pertStats, sname, uname, zname, printz)

% -------------------------------------------------------------------------
% prompt for shared limits, matching the old legacy behavior
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
    [xr, yr] = mc_fit_dialogs('ask_xy_limits', xr0, yr0);

    if ~isempty(xr) && ~isempty(yr)
        xlim(hxy(1),xr); xlim(xz,xr); xlim(yz,xr); xlim(hpert(1),xr);
        ylim(hxy(1),yr); ylim(xz,yr); ylim(yz,yr); ylim(hpert(1),yr);
    end
end

% -------------------------------------------------------------------------
% output naming and export behavior
file = strrep(file, "_perturbed.pts", "");

mc_fit_export_plot(gcf, file, outPDF);

function local_copy_plot_children(ax_from, ax_to, exclude_handles)
    if nargin < 3 || isempty(exclude_handles)
        exclude_handles = gobjects(0,1);
    end
    kids = allchild(ax_from);
    for k = numel(kids):-1:1
        hk = kids(k);
        if ~isgraphics(hk) || strcmp(get(hk,'Type'),'legend')
            continue
        end
        if any(hk == exclude_handles)
            continue
        end
        copyobj(hk, ax_to);
    end
end

function local_reset_perturbation_legend(ax, handles, stats, sname, uname, zname, printz)
    lgd = legend(ax);
    if isgraphics(lgd)
        delete(lgd)
    end

    hlist = gobjects(0,1);
    tlist = {};

    trial_handle = local_first_graphic([handles.trials_fit; handles.trials_nofit]);
    if isgraphics(trial_handle)
        hlist(end+1,1) = trial_handle; %#ok<AGROW>
        tlist{end+1} = local_make_perturbation_text(stats, sname, uname, printz); %#ok<AGROW>
    end

    if isgraphics(handles.overall)
        hlist(end+1,1) = handles.overall; %#ok<AGROW>
        tlist{end+1} = local_make_best_overall_text(stats, sname, uname, zname, printz); %#ok<AGROW>
    end

    if ~isempty(hlist)
        hleg = legend(ax, hlist, tlist);
        hleg.Location = 'northwest';
    end
end

function h = local_first_graphic(arr)
    h = gobjects(1,1);
    for i = 1:numel(arr)
        if isgraphics(arr(i))
            h = arr(i);
            return
        end
    end
end

function strg = local_make_perturbation_text(stats, sname, uname, printz)
    strg = "\bfPerturbations\rm" + newline ...
        + "\sigma_{" + sname(1) + "}\rm = \pm\rm" + num2str(stats.sigx,2) + uname(1) ...
        + newline + "\sigma_{" + sname(2) + "}\rm = \pm\rm" + num2str(stats.sigy,2) + uname(2);
    if printz
        strg = strg + newline + "\epsilon_{\rm" + sname(3) + "}\rm = \pm\rm" + num2str(stats.epsz,3);
    end
end

function strg = local_make_best_overall_text(stats, sname, uname, zname, printz)
    strg = "\bfBest overall\rm" + newline ...
        + sname(1) + "\rm = \rm" + num2str(stats.overall_x,4) + uname(1) + newline ...
        + sname(2) + "\rm = \rm" + num2str(stats.overall_y,4) + uname(2);
    if printz
        strg = strg + newline + "\rm" + zname + "\rm = \rm" + num2str(stats.overall_z,3);
    end
end
