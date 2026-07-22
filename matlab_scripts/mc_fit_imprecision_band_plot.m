% mc_fit_imprecision_interval_plot.m
% Plot the perturbation and corresponding central mc_fit *.pts files in an
% x-ln(Misfit) or y-ln(Misfit) projection, and draw the perturbation-derived
% imprecision interval on both axes.
%
% The first selected file must be the *_perturbed.pts file. The script then
% reads the corresponding *_central.pts file automatically, as in mc_fit_plot.
%
% The imprecision interval is defined as
%
%     ln(BCM_misfit) <= ln(Misfit) <= ln(BCM_misfit) + epsilon_ln(Misfit)
%
% where epsilon_ln(Misfit) is estimated from the perturbation analysis. In the
% central-model panel, central models outside this interval are also plotted,
% but with open symbols.

clear; close all; clc

% -------------------------------------------------------------------------
% switches: keep these aligned with mc_fit_plot except where noted
print   = true;
outPDF  = false;
plotMPP = false;
plotcov = false;   % required for this x/y-ln(Misfit) projection figure
ploterr = true;    % same error-bar legend logic as mc_fit_plot
plotbst = true;
auto_filled_symbols_for_nofit = true;
filled_symbols_for_nofit = false;
printz  = false;
maksub  = true;
fitonly_prompt = false;
manual_filter  = false;
filter         = true;    % central statistics/legend use mc_fit_plot filtering
% -------------------------------------------------------------------------
% Projection: 'xz' gives x-ln(Misfit); 'yz' gives y-ln(Misfit).
projection = 'xz';
ask_projection = true;
% -------------------------------------------------------------------------
% Unit conversions
convert_to_MPa  = false;
convert_to_kbar = true;
convert_to_GPa  = false;
convert_to_C    = false;
% -------------------------------------------------------------------------
% user parameters
sigma_level    = 1;
target_rse     = 0.20;
cap_fraction   = 0.10;
legend_coord_digits = 4;
legend_sigma_digits = 2;
outlier_filter = true;
outlier_filter_iterations = 10;
outlier_filter_convergence = 0.01;
outlier_filter_z = true;
outlier_filter_delta_z = false;
outlier_filter_delta_z_max = 3.91;
interval_label_factor = 1.1; %horizontal position relative to the BCM 
% -------------------------------------------------------------------------
% Interval rendering
interval_color = [1.0 0.85 0.0];
interval_alpha = 0.10;
arrow_line_width = get(groot,'DefaultLineLineWidth');
% -------------------------------------------------------------------------
% hardening
assert(exist('sigma_level','var') == 1, 'sigma_level not defined')
assert(exist('target_rse','var') == 1, 'target_rse not defined')
assert(exist('cap_fraction','var') == 1, 'cap_fraction not defined')
assert(exist('outlier_filter','var') == 1, 'outlier_filter not defined')
assert(exist('outlier_filter_iterations','var') == 1, 'outlier_filter_iterations not defined')
assert(exist('outlier_filter_convergence','var') == 1, 'outlier_filter_convergence not defined')
assert(exist('outlier_filter_z','var') == 1, 'outlier_filter_z not defined')
assert(exist('outlier_filter_delta_z','var') == 1, 'outlier_filter_delta_z not defined')
assert(exist('outlier_filter_delta_z_max','var') == 1, 'outlier_filter_delta_z_max not defined')

% -------------------------------------------------------------------------
% initialization
first = true;
function_to_set_graphics_defaults(print)
figure(1); clf
set(gcf,'Color','w')

if ask_projection
    choice = questdlg('Choose imprecision projection:', ...
        'MC_fit imprecision interval projection', ...
        'x-ln(Misfit)','y-ln(Misfit)','x-ln(Misfit)');
    if strcmp(choice,'y-ln(Misfit)')
        projection = 'yz';
    elseif strcmp(choice,'x-ln(Misfit)')
        projection = 'xz';
    elseif isempty(choice)
        error('No projection selected.')
    end
end

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
[xp,yp,zp,symp,fitp,xname,yname,zname,sname,uname,~,xvar,yvar,zvar,file,path,model,filled_symbols_for_nofit_auto,~] = ...
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

% Optional fit-based filtering, copied from mc_fit_plot.
filter_by_fit = false;
fitonly_prompt_central = false;
pert_sym = 2 + double(zname == "ln(Bayes)") * 3;
idx_pert_trials = find(symp == pert_sym & ~isnan(xp) & ~isnan(yp) & ~isnan(zp));
if ~isempty(idx_pert_trials)
    nfit = sum(fitp(idx_pert_trials) == 1);
    ntrial = numel(idx_pert_trials);
    if nfit > 0 && nfit < ntrial
        filter_by_fit = mc_fit_dialogs('ask_fit_filter_single', nfit, ntrial, 'perturbation');
        if filter_by_fit
            idx = idx_pert_trials(fitp(idx_pert_trials) ~= 1);
            xp(idx) = NaN; yp(idx) = NaN; zp(idx) = NaN;
        else
            fitonly_prompt_central = true;
        end
    end
end

[plot1, pertStats, ~, pcovar, limits] = function_to_do_mc_fit_plots( ...
    xp, yp, zp, symp, fitp, xname, yname, zname, sname, uname, model, filter, ...
    pcovar, limits, opts, {{1,2,1}});

pert_cache = struct('x', xp, 'y', yp, 'z', zp, 'symb', symp, 'fit', fitp, ...
    'xname', xname, 'yname', yname, 'zname', zname, ...
    'sname', sname, 'uname', uname, 'model', model);

% -------------------------------------------------------------------------
% switch to central file; read all central models, then apply the imprecision
% interval only to the arrays used for central statistics and the legend.
file = strrep(file, "perturbed", "central");
first = false;

[xc_all,yc_all,zc_all,symc_all,fitc_all,xname,yname,zname,sname,uname,~,xvar,yvar,zvar,file,path,model,filled_symbols_for_nofit_auto,~] = ...
    function_to_get_mc_fit_plot_file( ...
    '*_central.pts', xvar, yvar, zvar, NaN, ...
    fitonly_prompt_central, manual_filter, false, xname, yname, zname, sname, file, path, first, ...
    convert_to_MPa, convert_to_kbar, convert_to_GPa, convert_to_C, filter_by_fit, ...
    false, outlier_filter_iterations, outlier_filter_convergence, ...
    outlier_filter_z, outlier_filter_delta_z, outlier_filter_delta_z_max);

if auto_filled_symbols_for_nofit
    opts.filled_symbols_for_nofit = logical(filled_symbols_for_nofit_auto);
end
opts_central = opts;
opts_central.imprecision_filter_info = struct('applied', true, ...
    'pre_count', NaN, 'post_count', NaN, 'removed_count', NaN);

idx_bcm = find(symc_all == 1 & isfinite(xc_all) & isfinite(yc_all) & isfinite(zc_all), 1, 'first');
if isempty(idx_bcm)
    error('No best central model could be identified in the central file.')
end
if ~isfield(pertStats,'epsz') || ~isfinite(pertStats.epsz)
    error('Could not determine epsilon_ln(Misfit) from the perturbation data.')
end

interval_z0 = zc_all(idx_bcm);
interval_eps = pertStats.epsz;
opts_central.central_epsz_override = interval_eps;
interval_z1 = interval_z0 + interval_eps;

% Keep the best central model and central trials inside the interval for the
% filtered central analysis. Preserve all original fit flags and symbol colors.
bayes = double(zname == "ln(Bayes)") * 3;
cent_sym = 3 + bayes;
idx_cent_all = find(symc_all == cent_sym & isfinite(xc_all) & isfinite(yc_all) & isfinite(zc_all));
inside_interval = false(size(zc_all));
inside_interval(idx_bcm) = true;
inside_interval(idx_cent_all) = zc_all(idx_cent_all) >= interval_z0 & zc_all(idx_cent_all) <= interval_z1;

xc_filt = xc_all; yc_filt = yc_all; zc_filt = zc_all;
idx_hide = idx_cent_all(~inside_interval(idx_cent_all));
xc_filt(idx_hide) = NaN; yc_filt(idx_hide) = NaN; zc_filt(idx_hide) = NaN;

opts_central.imprecision_filter_info.pre_count = numel(idx_cent_all) + 1;
opts_central.imprecision_filter_info.post_count = sum(inside_interval([idx_bcm; idx_cent_all(:)]));
opts_central.imprecision_filter_info.removed_count = numel(idx_hide);

[plot2, centralStats, centralHandles, pcovar, limits] = function_to_do_mc_fit_plots( ...
    xc_filt, yc_filt, zc_filt, symc_all, fitc_all, xname, yname, zname, sname, uname, model, true, ...
    pcovar, limits, opts_central, {{1,2,2}});

% Overlay central models outside the interval, without adding legend entries.
% Symbol fill is governed only by the fit-all-data rule; membership in the
% imprecision interval is not encoded by marker filling.
plot_central_outside_interval_fit_symbols(plot2(1), xc_all, yc_all, zc_all, fitc_all, symc_all, idx_hide, zname, opts_central.filled_symbols_for_nofit)

% Redraw perturbation subplot after the central analysis is finalized. Force a
% filled BCM marker while keeping the mc_fit_plot legend text and symbol colors.
opts_redraw = opts_pert;
fit_for_filled_bcm = double(~opts_redraw.filled_symbols_for_nofit);
opts_redraw.best_override = struct('x', centralStats.best_x, ...
    'y', centralStats.best_y, 'z', centralStats.best_z, ...
    'fit', fit_for_filled_bcm);
if ~isempty(plot1) && isgraphics(plot1(1))
    cla(plot1(1))
end
[plot1, pertStats, ~, pcovar, limits] = function_to_do_mc_fit_plots( ...
    pert_cache.x, pert_cache.y, pert_cache.z, pert_cache.symb, pert_cache.fit, ...
    pert_cache.xname, pert_cache.yname, pert_cache.zname, ...
    pert_cache.sname, pert_cache.uname, pert_cache.model, filter, ...
    pcovar, limits, opts_redraw, {{1,2,1}});

% Overlay perturbation-derived data error bars onto the central-model plot,
% centered on the best central model, as in mc_fit_plot.
if ploterr && ~isempty(plot2) && isgraphics(plot2(1))
    hover = function_errorbar3_caps( ...
        plot2(1), centralStats.best_x, centralStats.best_y, centralStats.best_z, ...
        sigma_level * pertStats.sigx, sigma_level * pertStats.sigy, interval_eps, ...
        'Color', 'red', 'LineWidth', get(groot,'DefaultLineLineWidth'), ...
        'cap_fraction', cap_fraction);

    function_stack_errorbar_groups_xy(centralHandles.imprecision, hover, ...
        centralStats.sigx, centralStats.sigy, pertStats.sigx, pertStats.sigy);
end

% Set the requested projection before drawing the interval and label.
apply_projection(plot1(1), projection)
apply_projection(plot2(1), projection)

% Draw the interval on both plotted axes. The band spans the displayed axis by
% limit listeners, and the label/arrows are stored in data coordinates.
ax = [plot1(1), plot2(1)];
bcm_x = xc_all(idx_bcm);
bcm_y = yc_all(idx_bcm);
for iax = 1:2
    draw_imprecision_interval_on_axis(ax(iax), projection, interval_z0, interval_z1, ...
        interval_color, interval_alpha, ax(iax).LineWidth, interval_label_factor, ...
        pert_cache.x, pert_cache.y, pert_cache.z, xc_all, yc_all, zc_all, bcm_x, bcm_y);
end

% -------------------------------------------------------------------------
% final touches
% Use the coordinate range of both plotted data sets in the common-limits
% dialog.  The central plot includes central models outside the imprecision
% interval as open symbols, so include those coordinates as well.
valid_all = [isfinite(pert_cache.x(:)) & isfinite(pert_cache.y(:)) & isfinite(pert_cache.z(:)); ...
             isfinite(xc_all(:))       & isfinite(yc_all(:))       & isfinite(zc_all(:))];
x_all_limits = [pert_cache.x(:); xc_all(:)];
y_all_limits = [pert_cache.y(:); yc_all(:)];
z_all_limits = [pert_cache.z(:); zc_all(:)];
if any(valid_all)
    limits.minx = min(x_all_limits(valid_all)); limits.maxx = max(x_all_limits(valid_all));
    limits.miny = min(y_all_limits(valid_all)); limits.maxy = max(y_all_limits(valid_all));
    limits.minz = min(z_all_limits(valid_all)); limits.maxz = max(z_all_limits(valid_all));
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

file = strrep(file, "_central.pts", "_imprecision_interval");
mc_fit_export_plot(gcf, file, outPDF);

% =========================================================================
% Local helpers

function apply_projection(ax, projection)
    axes(ax) %#ok<LAXES>
    switch lower(projection)
        case 'xz'
            view(ax, 0, 0)     % x horizontal, z vertical
        case 'yz'
            view(ax, 90, 0)    % y horizontal, z vertical
        otherwise
            error('Unknown projection: %s', projection)
    end
    box(ax,'on')
end

function plot_central_outside_interval_fit_symbols(ax, x, y, z, fit, symb, idx, zname, filled_symbols_for_nofit)
    if isempty(idx) || ~isgraphics(ax,'axes')
        return
    end
    sc = ['b','r','b','r','b','r','b','r','b','r','b','r'];
    mt = ['o','o','s','s','v','v','d','d','s','s','v','v'];
    bayes = double(zname == "ln(Bayes)") * 3;
    cent_sym = 3 + bayes;
    idx = idx(symb(idx) == cent_sym & isfinite(x(idx)) & isfinite(y(idx)) & isfinite(z(idx)));
    if isempty(idx)
        return
    end

    hold(ax,'on')
    lw = 0.5*get(groot,'DefaultLineLineWidth');
    for fit_value = [1 0]
        jnd = idx(fit(idx) == fit_value);
        if isempty(jnd)
            continue
        end
        [fc, ec] = marker_fill_from_fit_rule(fit_value, sc(cent_sym), sc(cent_sym), filled_symbols_for_nofit);
        plot3(ax, x(jnd), y(jnd), z(jnd), mt(cent_sym), ...
            'LineStyle','none', 'MarkerFaceColor',fc, ...
            'MarkerEdgeColor', ec, ...
            'MarkerSize', 6, 'LineWidth', lw, ...
            'HandleVisibility','off');
    end
end

function [fc, ec] = marker_fill_from_fit_rule(fit_value, fill_color, edge_color, filled_symbols_for_nofit)
    % Match mc_fit_plot: marker fill encodes the fit-all-data flag only.
    % The imprecision interval is shown by the band and is not encoded by
    % open/filled central-model symbols.
    is_fit = (fit_value == 1);
    make_filled = xor(is_fit, logical(filled_symbols_for_nofit));
    if make_filled
        fc = fill_color;
    else
        fc = 'none';
    end
    ec = edge_color;
end

function draw_imprecision_interval_on_axis(ax, projection, z0, z1, face_color, face_alpha, lw, interval_label_factor, xp, yp, zp, xc, yc, zc, bcm_x, bcm_y)
    if ~isgraphics(ax,'axes') || ~all(isfinite([z0 z1]))
        return
    end
    hold(ax,'on')

    xl = xlim(ax); yl = ylim(ax);
    xspan = diff(xl); yspan = diff(yl);
    if xspan <= 0 || yspan <= 0
        return
    end

    switch lower(projection)
        case 'xz'
            yplane = yl(1) + 0.50*yspan;
            hband = patch(ax, [xl(1) xl(2) xl(2) xl(1)], [yplane yplane yplane yplane], ...
                [z0 z0 z1 z1], face_color, ...
                'FaceAlpha', face_alpha, 'EdgeColor', 'k', 'LineWidth', lw, ...
                'HandleVisibility','off', 'Clipping','on');
            try, uistack(hband,'bottom'); catch, end
            % horizontal coordinate, the 1st method is better if the label
            % is to be placed at the same position if common scales are used 
            hcoord = interval_label_factor * bcm_x;
            %hcoord = bcm_x + 0.1*diff(xlim(ax));
            if ~isfinite(hcoord)
                hcoord = choose_low_density_coordinate([xp(:); xc(:)], [zp(:); zc(:)], xl, [z0 z1]);
            end
            label_and_arrows = draw_interval_label_and_arrows(ax, hcoord, yplane, z0, z1, ...
                face_color, face_alpha, lw, 'xz');
            add_interval_limit_listener(ax, hband, label_and_arrows, projection, z0, z1)
        case 'yz'
            xplane = xl(1) + 0.50*xspan;
            hband = patch(ax, [xplane xplane xplane xplane], [yl(1) yl(2) yl(2) yl(1)], ...
                [z0 z0 z1 z1], face_color, ...
                'FaceAlpha', face_alpha, 'EdgeColor', 'k', 'LineWidth', lw, ...
                'HandleVisibility','off', 'Clipping','on');
            try, uistack(hband,'bottom'); catch, end
            % horizontal coordinate, the 1st method is better if the label
            % is to be placed at the same position if common scales are used 
            hcoord = interval_label_factor * bcm_y;
            %hcoord = bcm_y + 0.1*diff(ylim(ax));
            if ~isfinite(hcoord)
                hcoord = choose_low_density_coordinate([yp(:); yc(:)], [zp(:); zc(:)], yl, [z0 z1]);
            end
            label_and_arrows = draw_interval_label_and_arrows(ax, xplane, hcoord, z0, z1, ...
                face_color, face_alpha, lw, 'yz');
            add_interval_limit_listener(ax, hband, label_and_arrows, projection, z0, z1)
        otherwise
            error('Unknown projection: %s', projection)
    end
end

function add_interval_limit_listener(ax, hband, label_handles, projection, z0, z1)
    f = ancestor(ax,'figure');
    L = [addlistener(ax,'XLim','PostSet',@(~,~)update_interval_band(ax,hband,label_handles,projection,z0,z1)), ...
         addlistener(ax,'YLim','PostSet',@(~,~)update_interval_band(ax,hband,label_handles,projection,z0,z1)), ...
         addlistener(ax,'ZLim','PostSet',@(~,~)update_interval_band(ax,hband,label_handles,projection,z0,z1))];
    oldL = getappdata(f,'mc_fit_imprecision_interval_listeners');
    setappdata(f,'mc_fit_imprecision_interval_listeners',[oldL L]);
end

function update_interval_band(ax, hband, label_handles, projection, z0, z1)
    if ~isgraphics(ax,'axes') || ~isgraphics(hband)
        return
    end
    xl = xlim(ax); yl = ylim(ax);
    switch lower(projection)
        case 'xz'
            yplane = yl(1) + 0.50*diff(yl);
            set(hband, 'XData', [xl(1) xl(2) xl(2) xl(1)], ...
                'YData', [yplane yplane yplane yplane], 'ZData', [z0 z0 z1 z1]);
            move_label_y_or_x(label_handles, [], yplane)
        case 'yz'
            xplane = xl(1) + 0.50*diff(xl);
            set(hband, 'XData', [xplane xplane xplane xplane], ...
                'YData', [yl(1) yl(2) yl(2) yl(1)], 'ZData', [z0 z0 z1 z1]);
            move_label_y_or_x(label_handles, xplane, [])
    end
end

function move_label_y_or_x(hs, newx, newy)
    for k = 1:numel(hs)
        if ~isgraphics(hs(k)), continue, end

        if isa(hs(k),'matlab.graphics.primitive.Text')
            % Text objects use Position, not XData/YData/ZData.
            pos = get(hs(k),'Position');
            if ~isempty(newx), pos(1) = newx; end
            if ~isempty(newy), pos(2) = newy; end
            set(hs(k),'Position',pos);

        elseif isprop(hs(k),'XData') && isprop(hs(k),'YData')
            % Patch and line objects use XData/YData. Preserve the
            % number of vertices/points and only move the projection plane.
            if ~isempty(newx)
                xd = get(hs(k),'XData');
                set(hs(k),'XData',newx*ones(size(xd)));
            end
            if ~isempty(newy)
                yd = get(hs(k),'YData');
                set(hs(k),'YData',newy*ones(size(yd)));
            end
        end
    end
end

function hcoord = choose_low_density_coordinate(h, z, hlim, zlim_interval)
    valid = isfinite(h) & isfinite(z) & h >= hlim(1) & h <= hlim(2) & z >= zlim_interval(1) & z <= zlim_interval(2);
    nb = 12;
    edges = linspace(hlim(1), hlim(2), nb+1);
    if any(valid)
        counts = histcounts(h(valid), edges);
        if nb > 4
            search_bins = 3:(nb-2);
        else
            search_bins = 1:nb;
        end
        [~, k0] = min(counts(search_bins));
        k = search_bins(k0);
        hcoord = 0.5*(edges(k) + edges(k+1));
    else
        hcoord = hlim(1) + 0.55*diff(hlim);
    end
end

function hs = draw_interval_label_and_arrows(ax, xtext, ytext, z0, z1, face_color, face_alpha, lw, projection)
    zmid = 0.5*(z0 + z1);
    label = '\epsilon_{ln(Misfit)}';

    ht = text(ax, xtext, ytext, zmid, label, ...
        'Interpreter','tex', 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
        'BackgroundColor','none', 'EdgeColor','none', 'Margin', 3, ...
        'HandleVisibility','off', 'Clipping','on');
    drawnow limitrate

    [z_text_bottom, z_text_top, hmin, hmax] = estimate_text_extent(ht, z0, z1, xtext, ytext, projection);
    switch lower(projection)
        case 'xz'
            hbox = patch(ax, [hmin hmax hmax hmin], [ytext ytext ytext ytext], ...
                [z_text_bottom z_text_bottom z_text_top z_text_top], face_color, ...
                'FaceAlpha', face_alpha, 'EdgeColor','none', 'HandleVisibility','off', 'Clipping','on');
            try, uistack(hbox,'bottom'); catch, end
            try, uistack(ht,'top'); catch, end
        case 'yz'
            hbox = patch(ax, [xtext xtext xtext xtext], [hmin hmax hmax hmin], ...
                [z_text_bottom z_text_bottom z_text_top z_text_top], face_color, ...
                'FaceAlpha', face_alpha, 'EdgeColor','none', 'HandleVisibility','off', 'Clipping','on');
            try, uistack(hbox,'bottom'); catch, end
            try, uistack(ht,'top'); catch, end
    end

    % Place arrowheads slightly inside the imprecision interval so that
    % the heads remain visible when they point to the interval edges.
    epsz = z1 - z0;
    % arrow ends
    z_upper_arrow = z0 + 0.95*epsz;
    z_lower_arrow = z0 + 0.05*epsz;

    hup = plot3(ax, [xtext xtext], [ytext ytext], [z_text_top z_upper_arrow], '-k', ...
        'LineWidth', lw, 'Marker','^', 'MarkerIndices',2, 'MarkerFaceColor','k', ...
        'MarkerSize',5, 'HandleVisibility','off', 'Clipping','on');
    hdn = plot3(ax, [xtext xtext], [ytext ytext], [z_text_bottom z_lower_arrow], '-k', ...
        'LineWidth', lw, 'Marker','v', 'MarkerIndices',2, 'MarkerFaceColor','k', ...
        'MarkerSize',5, 'HandleVisibility','off', 'Clipping','on');
    hs = [hbox ht hup hdn];
end

function [zbot, ztop, hmin, hmax] = estimate_text_extent(ht, z0, z1, xtext, ytext, projection)
    zmid = 0.5*(z0 + z1);
    epsz = z1 - z0;
    % arrow starts
    zbot = zmid - 0.18*epsz;
    ztop = zmid + 0.12*epsz;

    switch lower(projection)
        case 'xz'
            hmin = xtext - 0.06*diff(xlim(ancestor(ht,'axes')));
            hmax = xtext + 0.06*diff(xlim(ancestor(ht,'axes')));
        case 'yz'
            hmin = ytext - 0.06*diff(ylim(ancestor(ht,'axes')));
            hmax = ytext + 0.06*diff(ylim(ancestor(ht,'axes')));
    end
    try
        ext = get(ht,'Extent');
        if numel(ext) >= 4 && all(isfinite(ext(1:4))) && ext(3) > 0 && ext(4) > 0
            switch lower(projection)
                case 'xz'
                    hmin = ext(1); hmax = ext(1) + ext(3);
                    zbot = ext(2); ztop = ext(2) + ext(4);
                case 'yz'
                    hmin = ext(1); hmax = ext(1) + ext(3);
                    zbot = ext(2); ztop = ext(2) + ext(4);
            end
            % Reject absurd extents sometimes returned before 3-D text has settled.
            if ztop <= zbot || abs(ztop-zbot) > 0.8*abs(epsz)
                zbot = zmid - 0.18*epsz;
                ztop = zmid + 0.12*epsz;
            end
        end
    catch
    end
end
