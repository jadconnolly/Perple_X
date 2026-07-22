% mc_fit_single_plot makes a single plot of a user selected pts file
%
% Updated to match the current function interfaces used by mc_fit_plot.m.
%
% JADC Aug 27, 2025
% Updated Apr 11, 2026

%                                       set print to true for settings
%                                       appropriate for printing, set 
%                                       print to false for settings 
%                                       appropriate for projection (ppt). 
print = true;
outPDF = false;
legend_coord_digits = 4;
legend_sigma_digits = 2;   % significant digits for legend uncertainties
function_to_set_graphics_defaults(print)
figure(1); clf

%                                       initialization
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
subspec = {};

first = true;                           % prompt for variable choices etc
fitonly_prompt = true;                  % prompt whether to only plot models
                                        % that fit all data
manual_filter = true;                   % manual filter
filter = false;                         % auto filter to within imprecision
outlier_filter = true;                    % filter perturbation outliers before plotting/statistics
outlier_filter_iterations = 10;
outlier_filter_convergence = 0.01;
outlier_filter_z = true;       % robust one-sided z filter on ln(Misfit)/ln(Bayes); rejects only the bad tail
outlier_filter_delta_z = false;
outlier_filter_delta_z_max = 3.91;

%                                       unit conversions
convert_to_MPa  = false;
convert_to_kbar = true;
convert_to_GPa  = false;
convert_to_C    = true;

%                                       outlier filter options
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

%                                       open the *.pts file
[x, y, z, symb, fit, xname, yname, zname, sname, uname, ~, xvar, yvar, zvar, ...
    file, path, model, filled_symbols_for_nofit_auto] = ...
    function_to_get_mc_fit_plot_file (...
    '*.pts', xvar, yvar, zvar, epsz, fitonly_prompt, ...
    manual_filter, filter, xname, yname, zname, sname, file, path, first, ...
    convert_to_MPa, convert_to_kbar, convert_to_GPa, convert_to_C, false, ...
    outlier_filter, outlier_filter_iterations, outlier_filter_convergence, ...
    outlier_filter_z, outlier_filter_delta_z, outlier_filter_delta_z_max);

file_root = regexprep(char(file), '(_central|_perturbed)\.pts$', '');
title_prefix = string(regexprep(file_root, '\.[Pp][Tt][Ss]$', ''));

%                               on/off switches
plotMPP = false;                % plot MPP in central model
plotcov = true;                 % plot covariance ellipse
ploterr = true;                 % plot error bars on central model
plotbst = (model == "Perturbation"); % plot best overall model for perturbed input
printz  = true;
auto_filled_symbols_for_nofit = true;  % true: infer marker convention from the fit-all-data distribution
filled_symbols_for_nofit = false;        % used only when auto_filled_symbols_for_nofit is false; true fills non-fitting models
maksub  = false;                % make subplots
sigma_level  = 1;               % 1 -> 68%%-style scaling
target_rse   = 0.20;            % relative precision threshold for std estimate warnings
cap_fraction = 0.10;            % cap size relative to orthogonal error


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

[plot1, stats, ~, pcovar, limits] = function_to_do_mc_fit_plots( ...
    x, y, z, symb, fit, xname, yname, zname, sname, uname, model, filter, ...
    pcovar, limits, opts, subspec);



%                                        final touches:
if ~all(isfinite([limits.minx limits.maxx limits.miny limits.maxy limits.minz limits.maxz]))
    valid = isfinite(x) & isfinite(y) & isfinite(z);
    if any(valid)
        limits.minx = min(x(valid)); limits.maxx = max(x(valid));
        limits.miny = min(y(valid)); limits.maxy = max(y(valid));
        limits.minz = min(z(valid)); limits.maxz = max(z(valid));
    end
end

apply_limits = mc_fit_dialogs('ask_apply_xyz_limits', 'Modify X-Y-Z limits?');

if apply_limits
    xr0 = [limits.minx limits.maxx];
    yr0 = [limits.miny limits.maxy];
    zr0 = [limits.minz limits.maxz];
    [xr, yr, zr] = mc_fit_dialogs('ask_xyz_limits', xr0, yr0, zr0);

    if ~isempty(xr) && ~isempty(yr) && ~isempty(zr)
        xlim(plot1(1), xr);
        ylim(plot1(1), yr);
        zlim(plot1(1), zr);
    end
end

file = regexprep(char(file), '(_central|_perturbed)\.pts$', '');
mc_fit_export_plot(gcf, file, outPDF);


