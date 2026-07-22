% 
% Perform the mc_fit_plot analysis for multiple problems and
% superimpose their uncertainty results on x-y-z subplots.
%
% At startup, select all *_perturbed.pts files and enter their BCM tags.
% The x-y inversion parameters are selected independently for every problem;
% the z coordinate is identified from the 999.000 objective-function marker.
%
% The model populations, best overall models, and error bars are not plotted.
% Each BCM and its total ellipse use one of 12 distinctive colors; colors are
% recycled after the twelfth problem. Total ellipses are outlined only.
%
% The legend contains only the named BCMs.  For ln(Misfit), the displayed
% total uncertainty is defined here by quadrature of the perturbation-derived
% imprecision and the retained-central-model z scatter:
%
%   epsilon_z,total = sqrt(epsilon_z,data^2 + epsilon_z,aleatoric^2)
%
% JADC / ChatGPT, July 2026

clear; close all; clc

% -------------------------------------------------------------------------
% User options (kept parallel to mc_fit_plot where applicable)
print   = true;
outPDF  = false;
filter  = true;
fitonly_prompt = false;
manual_filter  = false;

auto_filled_symbols_for_nofit = true;
filled_symbols_for_nofit = false;

convert_to_MPa  = false;
convert_to_kbar = true;
convert_to_GPa  = false;
convert_to_C    = false;

sigma_level = 1;
target_rse = 0.20;
legend_coord_digits = 4;
legend_sigma_digits = 2;

outlier_filter = true;
outlier_filter_iterations = 10;
outlier_filter_convergence = 0.01;
outlier_filter_z = true;
outlier_filter_delta_z = false;
outlier_filter_delta_z_max = 3.91;

total_ellipse_alpha = 0.0;
outfile_base = 'mc_fit_permutation_plot';

% Standard component colors used by mc_fit_plot.
color_data = [1 0 0];
color_aleatoric = [0 0 1];

% Twelve distinctive BCM/total colors.  Recycled cyclically if necessary.
bcm_colors = [ ...
    0.0000 0.4470 0.7410; ... % blue
    0.8500 0.3250 0.0980; ... % orange
    0.9290 0.6940 0.1250; ... % yellow
    0.4940 0.1840 0.5560; ... % purple
    0.4660 0.6740 0.1880; ... % green
    0.3010 0.7450 0.9330; ... % cyan
    0.3500 0.3500 0.7000; ... % slate blue
    0.2500 0.2500 0.2500; ... % dark gray
    0.6000 0.4000 0.2000; ... % brown
    0.9000 0.2000 0.6000; ... % magenta
    0.2000 0.6500 0.5000; ... % teal
    0.5000 0.5000 0.0000];    % olive

% -------------------------------------------------------------------------
% Initialization
function_to_set_graphics_defaults(print)

% Analysis results are accumulated first, then rendered into one three-panel
% figure and three standalone figures.
fig = [];
ax = [];

xvar = [];
yvar = [];
zvar = [];
epsz = NaN;
xname = "";
yname = "";
zname = "";
sname = strings(1,3);
uname = strings(1,3);

problem_count = 0;
records = struct('name',{},'color',{},'best_x',{},'best_y',{},'best_z',{}, ...
    'best_fit',{},'fill_rule',{},'Cdata',{},'Caleatoric',{},'Ctotal',{}, ...
    'sig_data',{},'sig_aleatoric',{},'sig_total',{});

% Collect the complete input set before any analysis or plotting.  Each
% uigetfile call may browse to a different directory.  After at least one
% file has been selected, Cancel ends the collection step.
selected_files = {};
selected_paths = {};
start_path = pwd;
while true
    if isempty(selected_files)
        prompt = sprintf('Select perturbation file %d (press Cancel when done)',1);
    else
        last_loaded = fullfile(selected_paths{end}, selected_files{end});
        prompt = sprintf(['Select perturbation file %d (press Cancel when done)' ...
            '\nLast file loaded:\n%s'], numel(selected_files)+1, last_loaded);
    end

    [selected_file, selected_path, indx] = uigetfile( ...
        {'*_perturbed.pts','MC_fit perturbation files (*_perturbed.pts)'}, ...
        prompt, start_path);

    if indx == 0
        if isempty(selected_files)
            error('No mc_fit problems were selected.')
        end
        break
    end

    selected_files{end+1,1} = selected_file; %#ok<SAGROW>
    selected_paths{end+1,1} = selected_path; %#ok<SAGROW>
    start_path = selected_path;
end

ndefault = numel(selected_files);
prompts = cell(ndefault,1);
default_tags = cell(ndefault,1);
for i = 1:ndefault
    full_input_name = fullfile(selected_paths{i}, selected_files{i});
    prompts{i} = sprintf('BCM tag for problem %d (%s):', i, full_input_name);
    default_tags{i} = sprintf('Problem %d', i);
end
tags = inputdlg(prompts, 'BCM tags', 1, default_tags);
if isempty(tags)
    error('No BCM tags were entered.')
end
for i = 1:numel(tags)
    tags{i} = strtrim(tags{i});
    if isempty(tags{i})
        tags{i} = default_tags{i};
    end
end

for iproblem = 1:numel(selected_files)
    file = selected_files{iproblem};
    path = selected_paths{iproblem};
    first_problem = (iproblem == 1);

    % Select x and y independently for every problem. Resetting the column
    % indices prevents a selection from a problem with a different number of
    % inversion parameters from being applied to this file. The file reader
    % identifies z from the 999.000 objective-function marker.
    xvar = [];
    yvar = [];
    zvar = [];
    select_variables = true;
    try
        [xp,yp,zp,symp,fitp,xname_p,yname_p,zname_p,sname_p,uname_p,~, ...
            xvar,yvar,zvar,file,path,modelp,filled_auto_p,~] = ...
            function_to_get_mc_fit_plot_file( ...
            '*_perturbed.pts', xvar, yvar, zvar, epsz, fitonly_prompt, ...
            manual_filter, filter, xname, yname, zname, sname, file, path, select_variables, ...
            convert_to_MPa, convert_to_kbar, convert_to_GPa, convert_to_C, false, ...
            outlier_filter, outlier_filter_iterations, outlier_filter_convergence, ...
            outlier_filter_z, outlier_filter_delta_z, outlier_filter_delta_z_max);
    catch ME
        rethrow(ME)
    end

    if modelp ~= "Perturbation"
        error('The selected file must be an *_perturbed.pts file.')
    end

    if first_problem
        xname = xname_p; yname = yname_p; zname = zname_p;
        sname = sname_p; uname = uname_p;
    else
        assert(xname_p == xname && yname_p == yname && zname_p == zname, ...
            'All problems must use the same plotted variables.')
    end

    bcm_name = tags{iproblem};

    % Match mc_fit_plot fit-filter workflow.
    filter_by_fit = false;
    fitonly_prompt_central = false;
    pert_sym = 2 + double(zname_p == "ln(Bayes)") * 3;
    idx_pert_trials = find(symp == pert_sym & isfinite(xp) & isfinite(yp) & isfinite(zp));
    if ~isempty(idx_pert_trials)
        nfit = sum(fitp(idx_pert_trials) == 1);
        ntrial = numel(idx_pert_trials);
        if nfit > 0 && nfit < ntrial
            filter_by_fit = mc_fit_dialogs('ask_fit_filter_single', ...
                nfit, ntrial, 'perturbation');
            if filter_by_fit
                idx = idx_pert_trials(fitp(idx_pert_trials) ~= 1);
                xp(idx) = NaN; yp(idx) = NaN; zp(idx) = NaN;
            else
                fitonly_prompt_central = true;
            end
        end
    end

    fill_rule = logical(filled_symbols_for_nofit);
    if auto_filled_symbols_for_nofit
        fill_rule = logical(filled_auto_p);
    end

    % Run the shared analysis in an invisible temporary axes.  This preserves
    % mc_fit_plot statistics, filtering, warnings, and limiting-case behavior,
    % while the requested overlay is drawn separately below.
    tmpfig = figure('Visible','off');
    cleanup_tmp = onCleanup(@() close_if_graphics(tmpfig));
    tmpax = axes(tmpfig);
    limits_tmp = struct('minx',Inf,'maxx',-Inf,'miny',Inf,'maxy',-Inf,'minz',Inf,'maxz',-Inf);
    pcovar = [];
    opts = struct( ...
        'plotMPP', false, 'plotcov', true, 'ploterr', false, ...
        'plotbst', false, 'printz', false, 'print', print, 'maksub', false, ...
        'target_axes', tmpax, ...
        'filled_symbols_for_nofit', fill_rule, ...
        'sigma_level', sigma_level, ...
        'legend_coord_digits', legend_coord_digits, ...
        'legend_sigma_digits', legend_sigma_digits, ...
        'cap_fraction', 0.10, 'target_rse', target_rse, ...
        'outlier_filter', logical(outlier_filter), ...
        'outlier_filter_iterations', outlier_filter_iterations, ...
        'outlier_filter_convergence', outlier_filter_convergence, ...
        'outlier_filter_z', logical(outlier_filter_z), ...
        'outlier_filter_delta_z', logical(outlier_filter_delta_z), ...
        'outlier_filter_delta_z_max', outlier_filter_delta_z_max);

    [~, pertStats, ~, pcovar, limits_tmp] = function_to_do_mc_fit_plots( ...
        xp, yp, zp, symp, fitp, xname_p, yname_p, zname_p, sname_p, uname_p, ...
        modelp, filter, pcovar, limits_tmp, opts, {});

    % Preserve the perturbation-derived data covariance.  The shared central
    % call returns the central covariance in the same output argument.
    Cdata = pcovar;

    % Read matching central file and apply perturbation-derived imprecision.
    central_file = strrep(file, 'perturbed', 'central');
    [xc,yc,zc,symc,fitc,xname_c,yname_c,zname_c,sname_c,uname_c,~, ...
        xvar,yvar,zvar,central_file,path,modelc,filled_auto_c,imprecision_filter_info] = ...
        function_to_get_mc_fit_plot_file( ...
        '*_central.pts', xvar, yvar, zvar, pertStats.filter_half_width_z, ...
        fitonly_prompt_central, manual_filter, filter, xname, yname, zname, ...
        sname, central_file, path, false, ...
        convert_to_MPa, convert_to_kbar, convert_to_GPa, convert_to_C, filter_by_fit, ...
        outlier_filter, outlier_filter_iterations, outlier_filter_convergence, ...
        outlier_filter_z, outlier_filter_delta_z, outlier_filter_delta_z_max);

    if modelc ~= "Central"
        error('The matching file must be an *_central.pts file.')
    end

    if auto_filled_symbols_for_nofit
        fill_rule = logical(filled_auto_c);
    end
    opts.filled_symbols_for_nofit = fill_rule;
    opts.imprecision_filter_info = imprecision_filter_info;
    cla(tmpax)

    [~, centralStats, ~, pcovar] = function_to_do_mc_fit_plots( ...
        xc, yc, zc, symc, fitc, xname_c, yname_c, zname_c, sname_c, uname_c, ...
        modelc, filter, pcovar, limits_tmp, opts, {});

    close_if_graphics(tmpfig)
    clear cleanup_tmp

    problem_count = problem_count + 1;
    bcm_color = bcm_colors(mod(problem_count-1,size(bcm_colors,1))+1,:);

    % Covariance matrices after all requested filtering.
    idx_cent = centralStats.idx_sigma_cov;
    if centralStats.aleatoric_covariance_available && numel(idx_cent) >= 3
        Caleatoric = cov([xc(idx_cent), yc(idx_cent)], 'omitrows');
    else
        Caleatoric = zeros(2);
    end
    Ctotal = Cdata + Caleatoric;

    % Store the complete plotting record. Rendering is deferred until all
    % problems have been analyzed so the same results can be drawn into the
    % combined and standalone figures.
    center_xy = [centralStats.best_x centralStats.best_y];

    sig_data = [ ...
        sigma_level * sqrt(max(Cdata(1,1),0)), ...
        sigma_level * sqrt(max(Cdata(2,2),0)), ...
        sigma_level * max(pertStats.epsz,0)];
    sig_aleatoric = [ ...
        sigma_level * sqrt(max(Caleatoric(1,1),0)), ...
        sigma_level * sqrt(max(Caleatoric(2,2),0)), ...
        sigma_level * max(centralStats.epsz,0)];
    sig_total = sqrt(sig_data.^2 + sig_aleatoric.^2);

    legend_name = bcm_name;
    if centralStats.best_fit == 1
        legend_name = [legend_name '*'];
    end

    records(problem_count).name = legend_name;
    records(problem_count).color = bcm_color;
    records(problem_count).best_x = center_xy(1);
    records(problem_count).best_y = center_xy(2);
    records(problem_count).best_z = centralStats.best_z;
    records(problem_count).best_fit = centralStats.best_fit;
    records(problem_count).fill_rule = fill_rule;
    records(problem_count).Cdata = Cdata;
    records(problem_count).Caleatoric = Caleatoric;
    records(problem_count).Ctotal = Ctotal;
    records(problem_count).sig_data = sig_data;
    records(problem_count).sig_aleatoric = sig_aleatoric;
    records(problem_count).sig_total = sig_total;
end


% -------------------------------------------------------------------------
% Common axis limits for all panels/figures.
allx = reshape([records.best_x],[],1);
ally = reshape([records.best_y],[],1);
allz = reshape([records.best_z],[],1);
for k = 1:numel(records)
    xy = covariance_ellipse_xy(records(k).Ctotal, ...
        [records(k).best_x records(k).best_y], sigma_level);
    if ~isempty(xy)
        allx = [allx; xy(:,1)]; %#ok<AGROW>
        ally = [ally; xy(:,2)]; %#ok<AGROW>
    end
end
common_limits = [min(allx) max(allx) min(ally) max(ally) min(allz) max(allz)];

% -------------------------------------------------------------------------
% Render one three-panel figure and three standalone figures.
fig_combined = figure(1); clf(fig_combined)
set(fig_combined,'Color','w')
tl = tiledlayout(fig_combined,1,3,'TileSpacing','compact','Padding','compact');
ax_total = nexttile(tl,1);
ax_data = nexttile(tl,2);
ax_aleatoric = nexttile(tl,3);

render_uncertainty_panel(ax_total, records, 'total', 'Total uncertainty', ...
    xname, yname, sname, uname, zname, sigma_level, total_ellipse_alpha, ...
    legend_coord_digits, legend_sigma_digits, common_limits);
render_uncertainty_panel(ax_data, records, 'data', 'Data uncertainty', ...
    xname, yname, sname, uname, zname, sigma_level, total_ellipse_alpha, ...
    legend_coord_digits, legend_sigma_digits, common_limits);
render_uncertainty_panel(ax_aleatoric, records, 'aleatoric', 'Aleatoric uncertainty', ...
    xname, yname, sname, uname, zname, sigma_level, total_ellipse_alpha, ...
    legend_coord_digits, legend_sigma_digits, common_limits);

fig_total = figure('Visible','off'); clf(fig_total); set(fig_total,'Color','w')
render_uncertainty_panel(axes(fig_total), records, 'total', 'Total uncertainty', ...
    xname, yname, sname, uname, zname, sigma_level, total_ellipse_alpha, ...
    legend_coord_digits, legend_sigma_digits, common_limits);

fig_data = figure('Visible','off'); clf(fig_data); set(fig_data,'Color','w')
render_uncertainty_panel(axes(fig_data), records, 'data', 'Data uncertainty', ...
    xname, yname, sname, uname, zname, sigma_level, total_ellipse_alpha, ...
    legend_coord_digits, legend_sigma_digits, common_limits);

fig_aleatoric = figure('Visible','off'); clf(fig_aleatoric); set(fig_aleatoric,'Color','w')
render_uncertainty_panel(axes(fig_aleatoric), records, 'aleatoric', 'Aleatoric uncertainty', ...
    xname, yname, sname, uname, zname, sigma_level, total_ellipse_alpha, ...
    legend_coord_digits, legend_sigma_digits, common_limits);

mc_fit_export_plot(fig_combined, [outfile_base '_combined'], outPDF)
mc_fit_export_plot(fig_total, [outfile_base '_total'], outPDF)
mc_fit_export_plot(fig_data, [outfile_base '_data'], outPDF)
mc_fit_export_plot(fig_aleatoric, [outfile_base '_aleatoric'], outPDF)


% =========================================================================
function render_uncertainty_panel(ax, records, kind, panel_title, ...
        xname, yname, sname, uname, zname, sigma_level, alpha_value, ...
        coord_digits, sigma_digits, common_limits)
    hold(ax,'on')
    box(ax,'on')
    axis(ax,'square')
    view(ax,0,90)
    xlabel(ax,xname)
    ylabel(ax,yname)
    title(ax,panel_title)

    n = numel(records);
    handles = gobjects(n,1);
    labels = cell(n,1);
    lw = get(groot,'DefaultLineLineWidth');

    for i = 1:n
        r = records(i);
        switch kind
            case 'total'
                C = r.Ctotal;
                sig = r.sig_total;
                use_patch = false;
            case 'data'
                C = r.Cdata;
                sig = r.sig_data;
                use_patch = false;
            case 'aleatoric'
                C = r.Caleatoric;
                sig = r.sig_aleatoric;
                use_patch = false;
            otherwise
                error('Unknown uncertainty panel type: %s',kind)
        end

        center_xy = [r.best_x r.best_y];
        if use_patch
            plot_covariance_patch(ax,C,center_xy,r.best_z,sigma_level, ...
                r.color,alpha_value,lw);
        else
            plot_covariance_outline(ax,C,center_xy,r.best_z,sigma_level, ...
                r.color,lw);
        end

        [fc,ec] = marker_fill_from_fit_local(r.best_fit,r.color,'k',r.fill_rule);
        handles(i) = plot3(ax,r.best_x,r.best_y,r.best_z,'o', ...
            'MarkerFaceColor',fc,'MarkerEdgeColor',ec, ...
            'MarkerSize',25,'LineWidth',lw);

        stats = struct('best_x',r.best_x,'best_y',r.best_y,'best_z',r.best_z);
        labels{i} = make_overlay_legend_text(r.name,stats,sig(1),sig(2),sig(3), ...
            sname,uname,zname,coord_digits,sigma_digits);
    end

    legend(ax,handles,labels,'Location','eastoutside','Interpreter','tex')
    xlim(ax,common_limits(1:2));
    ylim(ax,common_limits(3:4));
    zlim(ax,common_limits(5:6));
end

% =========================================================================
function plot_covariance_outline(ax, C, center_xy, z0, sigma_level, color, lw)
    xy = covariance_ellipse_xy(C, center_xy, sigma_level);
    if isempty(xy)
        return
    end
    plot3(ax, xy(:,1), xy(:,2), z0*ones(size(xy,1),1), '-', ...
        'Color', color, 'LineWidth', lw, 'HandleVisibility','off');
end

function plot_covariance_patch(ax, C, center_xy, z0, sigma_level, color, alpha_value, lw)
    xy = covariance_ellipse_xy(C, center_xy, sigma_level);
    if isempty(xy)
        return
    end
    patch(ax, xy(:,1), xy(:,2), z0*ones(size(xy,1),1), color, ...
        'FaceAlpha', alpha_value, 'EdgeColor', color, 'LineWidth', lw, ...
        'HandleVisibility','off');
end

function xy = covariance_ellipse_xy(C, center_xy, sigma_level)
    xy = [];
    if isempty(C) || any(size(C) ~= [2 2]) || any(~isfinite(C(:)))
        return
    end
    C = (C + C.')/2;
    [V,D] = eig(C);
    d = max(real(diag(D)),0);
    if all(d == 0)
        return
    end

    p = erf(sigma_level/sqrt(2));
    scale = sqrt(-2*log(1-p));
    theta = linspace(0,2*pi,181);
    circle = [cos(theta); sin(theta)];
    transform = V * diag(scale*sqrt(d));
    xy = (transform*circle).';
    xy(:,1) = xy(:,1) + center_xy(1);
    xy(:,2) = xy(:,2) + center_xy(2);
end

function [fc, ec] = marker_fill_from_fit_local(fit_value, fill_color, edge_color, filled_symbols_for_nofit)
    is_fit = (fit_value == 1);
    make_filled = xor(is_fit, filled_symbols_for_nofit);
    if make_filled
        fc = fill_color;
    else
        fc = 'none';
    end
    ec = edge_color;
end

function txt = make_overlay_legend_text(name, stats, sigx, sigy, sigz, ...
        sname, uname, zname, coord_digits, sigma_digits)
    xlab = legend_var_label(sname(1));
    ylab = legend_var_label(sname(2));
    zlab = plain_label(zname);
    xunit = plain_label(uname(1));
    yunit = plain_label(uname(2));

    txt = ['\bf' escape_tex(name) '\rm' newline ...
        '\rm' xlab ' = ' value_pm_unit(stats.best_x, sigx, xunit, coord_digits, sigma_digits) newline ...
        '\rm' ylab ' = ' value_pm_unit(stats.best_y, sigy, yunit, coord_digits, sigma_digits) newline ...
        '\rm' zlab ' = ' value_pm_unit(stats.best_z, sigz, '', coord_digits, sigma_digits)];
end

function out = value_pm_unit(value, sigma, unit, coord_digits, sigma_digits)
    out = [sigfig_local(value,coord_digits) ' \pm ' sigfig_local(sigma,sigma_digits)];
    if ~isempty(unit)
        out = [out ' ' unit];
    end
end

function out = sigfig_local(value,n)
    if isempty(value) || ~isfinite(value)
        out = num2str(value);
    else
        out = sprintf('%.*g',n,value);
    end
end

function out = plain_label(in)
    out = char(string(in));
    out = regexprep(out, '\\(it|rm|bf)\s*', '');
    out = strtrim(out);
end

function out = legend_var_label(in)
    base = plain_label(in);
    if contains(base, '_')
        parts = split(string(base), '_');
        head = char(parts(1));
        tail = char(join(parts(2:end), '_'));
        out = ['\it' head '_{\rm' tail '}'];
    else
        out = ['\it' base '\rm'];
    end
end

function out = escape_tex(in)
    out = char(string(in));
    out = strrep(out, '\', '\\');
    out = strrep(out, '_', '\_');
    out = strrep(out, '%', '\%');
end

function close_if_graphics(h)
    if ~isempty(h) && isgraphics(h)
        close(h)
    end
end
