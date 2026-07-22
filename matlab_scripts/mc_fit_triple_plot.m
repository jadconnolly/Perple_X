% legacy_mc_fit_tri_plot_script.m
% Three-panel version of legacy_mc_fit_plot as a plain SCRIPT.
%
% Behavior:
% - lets the user choose the input file first
% - infers central vs perturbed from the selected file name / model tag
% - plots only:
%       (1) large XY panel
%       (2) small XZ panel
%       (3) small YZ panel
% - if the selected file is central:
%       * retains legacy subplot-1 behavior
%       * allows manual filtering
%       * shows MPP in subplot 1
% - if the selected file is perturbed:
%       * suppresses MPP
%       * keeps the "best central model" blurb in subplot 1
%
% This version avoids local functions for broader MATLAB compatibility.

clf

print = true;
outPDF = false;
legend_coord_digits = 4;
legend_sigma_digits = 2;   % false: png bitmap; true: vector pdf
function_to_set_graphics_defaults(print)
figure(1); clf

% -------------------------------------------------------------------------
% switches
plotMPP = true;
plotcov = true;
ploterr = false;
plotbst = true;
auto_filled_symbols_for_nofit = true;  % true: infer marker convention from the fit-all-data distribution
filled_symbols_for_nofit = false;        % used only when auto_filled_symbols_for_nofit is false; true fills non-fitting models
printz  = true;
maksub  = true;

%                                       unit conversions
convert_to_MPa  = false;
convert_to_kbar = true;
convert_to_GPa  = false;
convert_to_C    = false;

sigma_level    = 1;
cap_fraction   = 0.10;
target_rse     = 0.20;
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
filter = false;
manual_filter = false;
pcovar = [];
limits = struct('minx',Inf,'maxx',-Inf,'miny',Inf,'maxy',-Inf,'minz',Inf,'maxz',-Inf);

% -------------------------------------------------------------------------
% First pass: choose the file and the plotting variables.
type = '*_*.pts';

[x, y, z, symb, fit, xname, yname, zname, sname, uname, ~, xvar, yvar, zvar, ...
    file, path, model, filled_symbols_for_nofit_auto] = ...
    function_to_get_mc_fit_plot_file( ...
    type, xvar, yvar, zvar, epsz, fitonly_prompt, ...
    manual_filter, filter, xname, yname, zname, sname, file, path, first, ...
    convert_to_MPa, convert_to_kbar, convert_to_GPa, convert_to_C, false, ...
    outlier_filter, outlier_filter_iterations, outlier_filter_convergence, ...
    outlier_filter_z, outlier_filter_delta_z, outlier_filter_delta_z_max);

% Infer mode from returned model first, with filename as fallback.
isCentral = false;
isPerturbed = false;

if isstring(model) || ischar(model)
    if strcmpi(string(model),"Central")
        isCentral = true;
    elseif strcmpi(string(model),"Perturbation")
        isPerturbed = true;
    end
end

if ~(isCentral || isPerturbed)
    fl = lower(char(file));
    isCentral   = contains(fl,'central');
    isPerturbed = contains(fl,'perturb');
end

if ~(isCentral || isPerturbed)
    warning(['Could not determine whether the selected file is central or ' ...
             'perturbed. Assuming central.'])
    isCentral = true;
    isPerturbed = false;
end

% For central input, re-read with manual filtering enabled so the legacy
% interactive filtering dialog is preserved.
if isCentral
    manual_filter = true;
    plotMPP = true;

    first = false;
    [x, y, z, symb, fit, xname, yname, zname, sname, uname, ~, xvar, yvar, zvar, ...
        file, path, model, filled_symbols_for_nofit_auto] = ...
        function_to_get_mc_fit_plot_file( ...
        type, xvar, yvar, zvar, epsz, fitonly_prompt, ...
        manual_filter, filter, xname, yname, zname, sname, file, path, first, ...
        convert_to_MPa, convert_to_kbar, convert_to_GPa, convert_to_C, false, ...
        outlier_filter, outlier_filter_iterations, outlier_filter_convergence, ...
    outlier_filter_z, outlier_filter_delta_z, outlier_filter_delta_z_max);
else
    manual_filter = false;
    plotMPP = false;
end

file_root = regexprep(char(file), '(_central|_perturbed)\.pts$', '');
title_prefix = string(file_root);

if auto_filled_symbols_for_nofit
    filled_symbols_for_nofit = logical(filled_symbols_for_nofit_auto);
else
    filled_symbols_for_nofit = logical(filled_symbols_for_nofit);
end

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
    'outlier_filter_delta_z_max', outlier_filter_delta_z_max, ...
    'title_prefix', title_prefix);

% -------------------------------------------------------------------------
% figure and layout: 2 x 3
%
%   [ XY  XY | XZ ]
%   [ XY  XY | YZ ]

subspec = {{2,3,[1,2,4,5]}};

[hxy, stats, handles, pcovar, limits] = function_to_do_mc_fit_plots( ...
    x, y, z, symb, fit, xname, yname, zname, sname, uname, model, filter, ...
    pcovar, limits, opts, subspec);

drawnow

% -------------------------------------------------------------------------
% Create the legacy XZ and YZ side panels by copying the XY plot while
% suppressing covariance ellipses, total covariance, MPP, and text labels.

xz = subplot(2,3,3);
yz = subplot(2,3,6);

exclude_handles = gobjects(0,1);
if isstruct(handles)
    if isfield(handles,'mpp') && ~isempty(handles.mpp)
        try
            exclude_handles = [exclude_handles; handles.mpp(:)];
        catch
        end
    end
    if isfield(handles,'cov') && ~isempty(handles.cov)
        try
            exclude_handles = [exclude_handles; handles.cov(:)];
        catch
        end
    end
    if isfield(handles,'cov_total') && ~isempty(handles.cov_total)
        try
            exclude_handles = [exclude_handles; handles.cov_total(:)];
        catch
        end
    end
end

kids = allchild(hxy(1));
for k = numel(kids):-1:1
    hk = kids(k);

    if ~isgraphics(hk)
        continue
    end

    hk_type = '';
    try
        hk_type = get(hk,'Type');
    catch
    end

    if strcmp(hk_type,'legend') || strcmp(hk_type,'text')
        continue
    end

    skip_this = false;
    if ~isempty(exclude_handles)
        skip_this = any(hk == exclude_handles);
    end

    if ~skip_this
        try
            copyobj(hk, xz);
        catch
        end
    end
end

kids = allchild(hxy(1));
for k = numel(kids):-1:1
    hk = kids(k);

    if ~isgraphics(hk)
        continue
    end

    hk_type = '';
    try
        hk_type = get(hk,'Type');
    catch
    end

    if strcmp(hk_type,'legend') || strcmp(hk_type,'text')
        continue
    end

    skip_this = false;
    if ~isempty(exclude_handles)
        skip_this = any(hk == exclude_handles);
    end

    if ~skip_this
        try
            copyobj(hk, yz);
        catch
        end
    end
end

hold(xz,'on')
hold(yz,'on')

xlim(xz,[limits.minx limits.maxx])
ylim(xz,[limits.miny limits.maxy])
zlim(xz,[limits.minz limits.maxz])

xlim(yz,[limits.minx limits.maxx])
ylim(yz,[limits.miny limits.maxy])
zlim(yz,[limits.minz limits.maxz])

view(xz,0,0)
view(yz,90,0)

xlabel(xz,xname)
ylabel(xz,yname)
zlabel(xz,zname)
axis(xz,'square')
box(xz,'on')

xlabel(yz,xname)
ylabel(yz,yname)
zlabel(yz,zname)
axis(yz,'square')
box(yz,'on')

if isPerturbed
    function_warn_if_sd_sample_too_small( ...
        'perturbation analysis', stats.n_sigma_points, target_rse, ...
        'Increase number_of_perturbations proportionally');
end

% -------------------------------------------------------------------------
% Prompt for XY limits on the XY panel only; keep full 3D limits on the
% rotated side views.

if ~all(isfinite([limits.minx limits.maxx limits.miny limits.maxy limits.minz limits.maxz]))
    valid = isfinite(x) & isfinite(y) & isfinite(z);
    if any(valid)
        limits.minx = min(x(valid)); limits.maxx = max(x(valid));
        limits.miny = min(y(valid)); limits.maxy = max(y(valid));
        limits.minz = min(z(valid)); limits.maxz = max(z(valid));
    end
end

apply_limits = mc_fit_dialogs('ask_apply_xyz_limits', 'Set the same X-Y-Z limits for the plots?');

if apply_limits
    xr0 = [limits.minx limits.maxx];
    yr0 = [limits.miny limits.maxy];
    [xr, yr] = mc_fit_dialogs('ask_xy_limits', xr0, yr0);
    if isempty(xr) || isempty(yr)
        return
    end

    xlim(hxy(1),xr)
    ylim(hxy(1),yr)

    xlim(xz,[limits.minx limits.maxx])
    ylim(xz,[limits.miny limits.maxy])
    zlim(xz,[limits.minz limits.maxz])

    xlim(yz,[limits.minx limits.maxx])
    ylim(yz,[limits.miny limits.maxy])
    zlim(yz,[limits.minz limits.maxz])
end


file = regexprep(char(file), '(_central|_perturbed)\.pts$', '');
mc_fit_export_plot(gcf, file, outPDF);
