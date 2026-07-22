% mc_fit_multi_plot.m
% Plot any number of mc_fit *_central.pts or *_perturbed.pts files as
% laterally adjacent native MATLAB panels.
%
% This script avoids copying saved .fig files. Each panel is drawn directly
% by function_to_do_mc_fit_plots into a tiledlayout axes, so line widths,
% fonts, titles, labels, and legends are generated natively in their final
% positions.

clear; close all; clc

% -------------------------------------------------------------------------
% User options, intentionally parallel to mc_fit_single_plot.m
print = true;
outPDF = false; %#ok<NASGU>  % retained for consistency; export is done below
legend_coord_digits = 4;
legend_sigma_digits = 2;
function_to_set_graphics_defaults(print)

% Unit conversions
convert_to_MPa  = false;
convert_to_kbar = true;
convert_to_GPa  = false;
convert_to_C    = false;

% Filtering options
fitonly_prompt = true;
manual_filter = true;
filter = false;                         % no automatic central filtering in this script
outlier_filter = true;
outlier_filter_iterations = 10;
outlier_filter_convergence = 0.01;
outlier_filter_z = true;       % robust one-sided z filter on ln(Misfit)/ln(Bayes); rejects only the bad tail
outlier_filter_delta_z = false;
outlier_filter_delta_z_max = 3.91;

% Plot options
plotMPP = false;
plotcov = false;
ploterr = true;
printz  = true;
maksub  = false;                        % target axes are supplied explicitly
sigma_level  = 1;
target_rse   = 0.20;
cap_fraction = 0.10;
plotbst_default = true;                 % only used for perturbed input panels
auto_filled_symbols_for_nofit = true;   % true: infer marker convention from the fit-all-data distribution
filled_symbols_for_nofit = false;       % used only when auto_filled_symbols_for_nofit is false; true fills non-fitting models

% Output options
outfile_base = 'combined_mc_fit_multi_plot';
export_resolution = 600;
link_xy_axes = false;                   % set true if all panels should share x/y limits

% -------------------------------------------------------------------------
% Initialization
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
limits0 = struct('minx',Inf,'maxx',-Inf,'miny',Inf,'maxy',-Inf,'minz',Inf,'maxz',-Inf);
subspec = {};

panel_data = struct([]);
first = true;

% -------------------------------------------------------------------------
% Query files. The first call uses function_to_get_mc_fit_plot_file so that
% the usual variable-selection dialog is used. Subsequent files are selected
% with uigetfile and use the same x/y/z variable columns.
while true

    if first
        [x, y, z, symb, fit, xname, yname, zname, sname, uname, ~, xvar, yvar, zvar, ...
            file, path, model, filled_symbols_for_nofit_auto] = function_to_get_mc_fit_plot_file( ...
            '*.pts', xvar, yvar, zvar, epsz, fitonly_prompt, ...
            manual_filter, filter, xname, yname, zname, sname, file, path, first, ...
            convert_to_MPa, convert_to_kbar, convert_to_GPa, convert_to_C, false, ...
            outlier_filter, outlier_filter_iterations, outlier_filter_convergence, ...
            outlier_filter_z, outlier_filter_delta_z, outlier_filter_delta_z_max);
        first = false;
    else
        [file2, path2, indx] = uigetfile('*.pts', ...
            sprintf('Select panel %d mc_fit *.pts file', numel(panel_data)+1), path);
        if isequal(file2,0) || indx == 0
            break
        end
        file = file2;
        path = path2;

        [x, y, z, symb, fit, xname, yname, zname, sname, uname, ~, xvar, yvar, zvar, ...
            file, path, model, filled_symbols_for_nofit_auto] = function_to_get_mc_fit_plot_file( ...
            '*.pts', xvar, yvar, zvar, epsz, fitonly_prompt, ...
            manual_filter, filter, xname, yname, zname, sname, file, path, false, ...
            convert_to_MPa, convert_to_kbar, convert_to_GPa, convert_to_C, false, ...
            outlier_filter, outlier_filter_iterations, outlier_filter_convergence, ...
            outlier_filter_z, outlier_filter_delta_z, outlier_filter_delta_z_max);
    end

    file_root = regexprep(char(file), '(_central|_perturbed)\.pts$', '');
    title_prefix = string(regexprep(file_root, '\.[Pp][Tt][Ss]$', ''));

    k = numel(panel_data) + 1;
    panel_data(k).x = x; %#ok<SAGROW>
    panel_data(k).y = y;
    panel_data(k).z = z;
    panel_data(k).symb = symb;
    panel_data(k).fit = fit;
    panel_data(k).xname = xname;
    panel_data(k).yname = yname;
    panel_data(k).zname = zname;
    panel_data(k).sname = sname;
    panel_data(k).uname = uname;
    panel_data(k).model = model;
    panel_data(k).filter = filter;
    panel_data(k).pcovar = pcovar;
    panel_data(k).limits = limits0;
    panel_data(k).title_prefix = title_prefix;
    panel_data(k).file = file;
    panel_data(k).path = path;
    if auto_filled_symbols_for_nofit
        panel_data(k).filled_symbols_for_nofit = logical(filled_symbols_for_nofit_auto);
    else
        panel_data(k).filled_symbols_for_nofit = logical(filled_symbols_for_nofit);
    end

    answer = questdlg('Add another *.pts file?', ...
        'Continue?', 'Yes','No','Yes');
    if ~strcmp(answer,'Yes')
        break
    end
end

nfig = numel(panel_data);
if nfig < 1
    error('No mc_fit *.pts files selected.')
end

% -------------------------------------------------------------------------
% Native multi-panel plot.  For exactly two panels, use the same ordinary
% subplot geometry as mc_fit_plot.m rather than tiledlayout.  This gives the
% two-panel multi-plot the same panel size, margins, and spacing as the
% standard perturbation/central mc_fit_plot output.  For three or more panels,
% keep tiledlayout because subplot geometry becomes too cramped.
if nfig == 2
    outfig = figure(1);
    clf(outfig)
    set(outfig,'Color','w')
    use_mc_fit_plot_subplots = true;
    t = []; %#ok<NASGU>
else
    outfig = figure('Color','w', ...
        'Units','normalized', ...
        'Position',[0.03 0.15 0.94 0.65]);

    t = tiledlayout(outfig, 1, nfig, ...
        'TileSpacing','compact', ...
        'Padding','compact');
    use_mc_fit_plot_subplots = false;
end

ax = gobjects(nfig,1);
all_limits = limits0;

for i = 1:nfig
    if use_mc_fit_plot_subplots
        ax(i) = subplot(1,2,i,'Parent',outfig);
    else
        ax(i) = nexttile(t);
    end

    pd = panel_data(i);

    opts = struct( ...
        'plotMPP', plotMPP, ...
        'plotcov', plotcov, ...
        'ploterr', ploterr, ...
        'plotbst', plotbst_default && pd.model == "Perturbation", ...
        'printz',  printz, ...
        'print',   print, ...
        'maksub',  maksub, ...
        'filled_symbols_for_nofit', logical(pd.filled_symbols_for_nofit), ...
        'sigma_level', sigma_level, ...
        'legend_coord_digits', legend_coord_digits, ...
        'legend_sigma_digits', legend_sigma_digits, ...
        'cap_fraction', cap_fraction, ...
        'target_rse', target_rse, ...
        'title_prefix', pd.title_prefix, ...
        'target_axes', ax(i));

    [~, ~, ~, pcovar_i, lim_i] = function_to_do_mc_fit_plots( ...
        pd.x, pd.y, pd.z, pd.symb, pd.fit, ...
        pd.xname, pd.yname, pd.zname, pd.sname, pd.uname, ...
        pd.model, pd.filter, pd.pcovar, pd.limits, opts, subspec);

    panel_data(i).pcovar = pcovar_i;
    panel_data(i).limits = lim_i;

    all_limits.minx = min(all_limits.minx, lim_i.minx);
    all_limits.maxx = max(all_limits.maxx, lim_i.maxx);
    all_limits.miny = min(all_limits.miny, lim_i.miny);
    all_limits.maxy = max(all_limits.maxy, lim_i.maxy);
    all_limits.minz = min(all_limits.minz, lim_i.minz);
    all_limits.maxz = max(all_limits.maxz, lim_i.maxz);
end

if link_xy_axes && nfig > 1
    linkaxes(ax,'xy')
end

% Optional global x-y limit modification.
apply_limits = mc_fit_dialogs('ask_apply_xyz_limits', 'Modify X-Y-Z limits?');
if apply_limits
    xr0 = [all_limits.minx all_limits.maxx];
    yr0 = [all_limits.miny all_limits.maxy];
    zr0 = [all_limits.minz all_limits.maxz];
    [xr, yr, ~] = mc_fit_dialogs('ask_xyz_limits', xr0, yr0, zr0);
    if ~isempty(xr) && ~isempty(yr)
        for i = 1:nfig
            xlim(ax(i), xr);
            ylim(ax(i), yr);
        end
    end
end

% -------------------------------------------------------------------------
% Save/export
outdir = panel_data(1).path;
outfile = fullfile(outdir, outfile_base);

savefig(outfig, [outfile '.fig']);

exportgraphics(outfig, [outfile '.pdf'], ...
    'BackgroundColor','white', ...
    'ContentType','vector');

exportgraphics(outfig, [outfile '.png'], ...
    'BackgroundColor','white', ...
    'Resolution',export_resolution);

fprintf('Saved:\n%s.fig\n%s.pdf\n%s.png\n', outfile, outfile, outfile);
