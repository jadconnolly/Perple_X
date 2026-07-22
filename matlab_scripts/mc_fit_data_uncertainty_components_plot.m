% mc_fit_data_uncertainty_components_plot plots thermodynamic and analytical
% perturbation components and their independent sum.
%
% The first user-selected *_perturbed.pts file is interpreted as the
% thermodynamic component of the Total data uncertainty.  The second is interpreted
% as the analytical component.  The covariance matrices are added assuming
% independence to obtain the total data-uncertainty covariance.
%
% JADC / ChatGPT, Jul 7, 2026

print = true;
outPDF = false;
legend_coord_digits = 4;
show_BCM_in_legend = false;    % include BCM entry in legend
legend_sigma_digits = 2;
function_to_set_graphics_defaults(print)
figure(1); clf
ax = gca;
hold(ax,'on')

% -------------------------------------------------------------------------
% User options
plot_points = true;
plotcov = true;
ploterr = true;
sigma_level = 1;
target_rse = 0.20;
cap_fraction = 0.10;

% Colors requested by user.
color_thermo = [0.4940 0.1840 0.5560];  % purple
color_analyt = [0.4660 0.6740 0.1880];  % green
color_total  = [1.0000 0.0000 0.0000];  % red
color_bcm    = '#EDB120';               % standard BCM color used by mc_fit plots

% Filtering/reading options, matching mc_fit_single_plot defaults.
fitonly_prompt = true;
manual_filter = true;
filter = false;
outlier_filter = true;
outlier_filter_iterations = 10;
outlier_filter_convergence = 0.01;
outlier_filter_z = true;
outlier_filter_delta_z = false;
outlier_filter_delta_z_max = 3.91;

% Unit conversions.
convert_to_MPa  = false;
convert_to_kbar = true;
convert_to_GPa  = false;
convert_to_C    = false;

% Initial values passed to the reader.
xvar = [];
yvar = [];
zvar = [];
epsz = NaN;
xname = "";
yname = "";
zname = "";
sname = strings(1,3);
uname = strings(1,3);
file = '';
path = '';

% -------------------------------------------------------------------------
% Select/read the thermodynamic perturbation file.  Variable choices are made
% here and reused for the analytical perturbation file.
[fileT, pathT, indx] = select_component_file('thermodynamic');
if indx == 0
    error('You did not choose the thermodynamic perturbation file, I quit!')
end

[xT, yT, zT, symT, fitT, xname, yname, zname, sname, uname, ~, xvar, yvar, zvar, ...
    fileT, pathT, modelT, filled_symbols_for_nofit_auto] = ...
    function_to_get_mc_fit_plot_file( ...
    '*.pts', xvar, yvar, zvar, epsz, fitonly_prompt, ...
    manual_filter, filter, xname, yname, zname, sname, fileT, pathT, true, ...
    convert_to_MPa, convert_to_kbar, convert_to_GPa, convert_to_C, false, ...
    outlier_filter, outlier_filter_iterations, outlier_filter_convergence, ...
    outlier_filter_z, outlier_filter_delta_z, outlier_filter_delta_z_max);

if modelT ~= "Perturbation"
    error('The first file must be a perturbation-analysis *_perturbed.pts file.')
end

% Select/read the analytical perturbation file, but reuse the first file's
% x/y/z variable choices so both covariance matrices are in the same space.
[fileA, pathA, indx] = select_component_file('analytical');
if indx == 0
    error('You did not choose the analytical perturbation file, I quit!')
end

[xA, yA, zA, symA, fitA, ~, ~, znameA, ~, ~, ~, ~, ~, ~, ...
    fileA, pathA, modelA, filled_symbols_for_nofit_auto_A] = ...
    function_to_get_mc_fit_plot_file( ...
    '*.pts', xvar, yvar, zvar, epsz, fitonly_prompt, ...
    manual_filter, filter, xname, yname, zname, sname, fileA, pathA, false, ...
    convert_to_MPa, convert_to_kbar, convert_to_GPa, convert_to_C, false, ...
    outlier_filter, outlier_filter_iterations, outlier_filter_convergence, ...
    outlier_filter_z, outlier_filter_delta_z, outlier_filter_delta_z_max);

if modelA ~= "Perturbation"
    error('The second file must be a perturbation-analysis *_perturbed.pts file.')
end
if znameA ~= zname
    warning('The two files use different objective-function labels: %s and %s.', zname, znameA)
end

filled_symbols_for_nofit = logical(filled_symbols_for_nofit_auto | filled_symbols_for_nofit_auto_A);
dl = get(groot,'DefaultLineLineWidth');

% Symbol codes for perturbation rows depend on whether the objective is
% ln(Misfit) or ln(Bayes).
bayes = double(zname == "ln(Bayes)") * 3;
pert_sym = 2 + bayes;

S_T = component_stats(xT, yT, zT, symT, fitT, pert_sym, sigma_level);
S_A = component_stats(xA, yA, zA, symA, fitA, pert_sym, sigma_level);

% Use the first file's BCM as the common center.  The two perturbation files
% should normally have the same BCM; warn if they do not.
center = S_T.best;
if any(abs([S_A.best.x - center.x, S_A.best.y - center.y, S_A.best.z - center.z]) > 100*eps(max(1,abs([center.x center.y center.z]))))
    warning(['The two perturbation files do not have identical BCM coordinates. ', ...
        'The total covariance and all error bars/ellipses are centered on the BCM from the thermodynamic file.'])
end

C_total = S_T.C + S_A.C;
sig_total = sigma_level * sqrt(max(diag(C_total),0));

% Sample-size warning, equivalent to the standard plotting routine.
n_req = ceil(1 + 1/(2*target_rse^2));
warn_small_sample(S_T.n, n_req, target_rse, 'thermodynamic')
warn_small_sample(S_A.n, n_req, target_rse, 'analytical')

% -------------------------------------------------------------------------
% Plot perturbation clouds, BCM, component error bars/ellipses, and total.
if plot_points
    plot_component_points(ax, xT, yT, zT, fitT, S_T.idx, color_thermo, filled_symbols_for_nofit, dl)
    plot_component_points(ax, xA, yA, zA, fitA, S_A.idx, color_analyt, filled_symbols_for_nofit, dl)
end

hBCM = plot_model_marker(ax, center.x, center.y, center.z, S_T.best.fit, color_bcm, 'o', 25, dl, filled_symbols_for_nofit);
% If the BCM marker is filled, keep the standard fill color but use a black edge
% so it remains distinct from the covariance/error-bar colors.
if isgraphics(hBCM)
    bcm_face = get(hBCM, 'MarkerFaceColor');
    if ~(ischar(bcm_face) && strcmpi(bcm_face, 'none'))
        set(hBCM, 'MarkerEdgeColor', 'k');
    end
end

hT = gobjects(1,1);
hA = gobjects(1,1);
hTot = gobjects(1,1);

if ploterr
    function_errorbar3_caps(ax, center.x, center.y, center.z, S_T.sig(1), S_T.sig(2), 0, ...
        'Color', color_thermo, 'LineWidth', dl, 'cap_fraction', cap_fraction);
    function_errorbar3_caps(ax, center.x, center.y, center.z, S_A.sig(1), S_A.sig(2), 0, ...
        'Color', color_analyt, 'LineWidth', dl, 'cap_fraction', cap_fraction);
    function_errorbar3_caps(ax, center.x, center.y, center.z, sig_total(1), sig_total(2), 0, ...
        'Color', color_total, 'LineWidth', dl, 'cap_fraction', cap_fraction);
end

if plotcov
    hT = plot_covariance_matrix_ellipse(ax, S_T.C, center.x, center.y, center.z, color_thermo, dl, sigma_level);
    hA = plot_covariance_matrix_ellipse(ax, S_A.C, center.x, center.y, center.z, color_analyt, dl, sigma_level);
    hTot = plot_covariance_matrix_ellipse(ax, C_total, center.x, center.y, center.z, color_total, dl, sigma_level);
end

format_axes_local(ax, xname, yname, zname)

% Legend: optional BCM entry + thermodynamic + analytical + total data uncertainty.
% No misfit/imprecision entries are included.
tlist_components = { ...
    make_component_text('thermodynamic', S_T.sig, sname, uname, legend_sigma_digits), ...
    make_component_text('analytical',    S_A.sig, sname, uname, legend_sigma_digits), ...
    make_component_text('data',          sig_total, sname, uname, legend_sigma_digits)};

if show_BCM_in_legend
    hlist = [hBCM hT hA hTot];
    tlist = [{make_best_model_text(center, sname, uname, zname, legend_coord_digits)}, tlist_components];
else
    hlist = [hT hA hTot];
    tlist = tlist_components;
end

valid_legend = isgraphics(hlist);
leg = legend(ax, hlist(valid_legend), tlist(valid_legend));
leg.Location = 'northwest';

title(ax, 'Data uncertainty components')
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
ax.ZAxis.Exponent = 0;

% Apply optional limits.
valid = isfinite([xT; xA]) & isfinite([yT; yA]) & isfinite([zT; zA]);
if any(valid)
    ask_limits = mc_fit_dialogs('ask_apply_xyz_limits', 'Modify X-Y-Z limits?');
    if ask_limits
        xr0 = [min([xT(valid(1:numel(xT))); xA(valid(numel(xT)+1:end))]) max([xT(valid(1:numel(xT))); xA(valid(numel(xT)+1:end))])];
        yr0 = [min([yT(valid(1:numel(yT))); yA(valid(numel(yT)+1:end))]) max([yT(valid(1:numel(yT))); yA(valid(numel(yT)+1:end))])];
        zr0 = [min([zT(valid(1:numel(zT))); zA(valid(numel(zT)+1:end))]) max([zT(valid(1:numel(zT))); zA(valid(numel(zT)+1:end))])];
        [xr, yr, zr] = mc_fit_dialogs('ask_xyz_limits', xr0, yr0, zr0);
        if ~isempty(xr) && ~isempty(yr) && ~isempty(zr)
            xlim(ax, xr); ylim(ax, yr); zlim(ax, zr);
        end
    end
end

% Maintain a square x-y presentation after any optional limit changes.
axis(ax, 'square')

out_root = regexprep(char(fileT), '(_central|_perturbed)\.pts$', '');
out_root = [out_root '_data_uncertainty_components'];
mc_fit_export_plot(gcf, out_root, outPDF);

% =========================================================================
% Local helper functions


function [file, path, indx] = select_component_file(component)
    component = lower(string(component));
    switch component
        case "thermodynamic"
            ttl = 'Step 1 of 2: Select the THERMODYNAMIC perturbation-analysis *.pts file';
        case "analytical"
            ttl = 'Step 2 of 2: Select the ANALYTICAL perturbation-analysis *.pts file';
        otherwise
            ttl = 'Select perturbation-analysis *.pts file';
    end
    [file, path, indx] = uigetfile({'*.pts','mc_fit perturbation files (*.pts)';'*.*','All Files (*.*)'}, ttl);
end

function S = component_stats(x, y, z, symb, fit, pert_sym, sigma_level)
    idx_best = find(symb == 1 & isfinite(x) & isfinite(y) & isfinite(z), 1, 'first');
    if isempty(idx_best)
        error('No best central model could be identified in a perturbation file.')
    end
    idx = find(symb == pert_sym & isfinite(x) & isfinite(y) & isfinite(z));
    if numel(idx) < 3
        error('At least three finite perturbation models are required to calculate a covariance ellipse.')
    end
    xy = [x(idx), y(idx)];
    S.C = cov(xy, 'omitrows');
    S.sig = sigma_level * sqrt(max(diag(S.C),0));
    S.n = numel(idx);
    S.idx = idx;
    S.best = struct('x', x(idx_best), 'y', y(idx_best), 'z', z(idx_best), 'fit', fit(idx_best));
end

function warn_small_sample(n, n_req, target_rse, label)
    if n < n_req
        warning(['Only %d %s perturbation point(s) are available for the variance estimates; ', ...
            'at least %d are required for the estimated standard deviation to be accurate to within %.0f%%.'], ...
            n, label, n_req, 100*target_rse)
    end
end

function h = plot_covariance_matrix_ellipse(ax, C, X0, Y0, Z0, col, lw, sigma_level)
    C = (C + C.')/2;
    [V,D] = eig(C);
    d = max(diag(D), 0);
    [largest_eigenval, idx] = max(d);
    largest_eigenvec = V(:,idx);
    smallest_eigenval = min(d);
    phi = atan2(largest_eigenvec(2), largest_eigenvec(1));
    if phi < 0, phi = phi + 2*pi; end
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
    h = plot3(ax, rell(:,1)+X0, rell(:,2)+Y0, Z0*ones(size(rell,1),1), '-', ...
        'Color', col, 'LineWidth', lw);
end

function plot_component_points(ax, x, y, z, fit, idx, col, filled_symbols_for_nofit, dl)
    sz = 6;
    lw = 0.5*dl;
    idx_fit = idx(fit(idx)==1);
    idx_nofit = idx(fit(idx)==0);
    if ~isempty(idx_fit)
        plot_model_marker(ax, x(idx_fit), y(idx_fit), z(idx_fit), 1, col, 'o', sz, lw, filled_symbols_for_nofit, col);
    end
    if ~isempty(idx_nofit)
        plot_model_marker(ax, x(idx_nofit), y(idx_nofit), z(idx_nofit), 0, col, 'o', sz, lw, filled_symbols_for_nofit, col);
    end
end

function h = plot_model_marker(ax, x, y, z, fit, ec, mt, sz, lw, filled_symbols_for_nofit, fc_override)
    if nargin < 11
        fc_override = [];
    end
    if all(fit == 1)
        if filled_symbols_for_nofit
            fc = 'none';
        else
            fc = fc_override;
            if isempty(fc), fc = ec; end
        end
    else
        if filled_symbols_for_nofit
            fc = fc_override;
            if isempty(fc), fc = ec; end
        else
            fc = 'none';
        end
    end
    h = plot3(ax, x, y, z, mt, 'MarkerSize', sz, 'LineWidth', lw, ...
        'MarkerEdgeColor', ec, 'MarkerFaceColor', fc, 'LineStyle', 'none');
end

function format_axes_local(ax, xname, yname, zname)
    grid(ax,'on')
    box(ax,'on')
    view(ax,2)
    axis(ax,'square')
    xlabel(ax, xname)
    ylabel(ax, yname)
    zlabel(ax, zname)
end

function strg = make_best_model_text(best, sname, uname, zname, nfmt)
    xlab = legend_var_label(sname(1));
    ylab = legend_var_label(sname(2));
    xunit = legend_plain_label(uname(1));
    yunit = legend_plain_label(uname(2));
    zlab = legend_plain_label(zname);
    strg = ['\bfBest central model\rm' newline ...
        '\rm' xlab ' = ' value_with_unit(best.x, xunit, nfmt) newline ...
        '\rm' ylab ' = ' value_with_unit(best.y, yunit, nfmt) newline ...
        '\rm' zlab ' = ' sigfig(best.z, nfmt)];
end

function strg = make_component_text(tag, sigxy, sname, uname, nfmt)
    xlab = legend_var_label(sname(1));
    ylab = legend_var_label(sname(2));
    xunit = legend_plain_label(uname(1));
    yunit = legend_plain_label(uname(2));
    switch string(tag)
        case "thermodynamic"
            head = '\bfThermodynamic data uncertainty\rm';
            subtag = 'thermo';
        case "analytical"
            head = '\bfAnalytical data uncertainty\rm';
            subtag = 'analytical';
        otherwise
            head = '\bfTotal data uncertainty\rm';
            subtag = 'data';
    end
    strg = [head newline ...
        '\rm\sigma_{' component_subscript(xlab, subtag) '} = ' value_with_unit(sigxy(1), xunit, nfmt) newline ...
        '\rm\sigma_{' component_subscript(ylab, subtag) '} = ' value_with_unit(sigxy(2), yunit, nfmt)];
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

function out = legend_plain_label(in)
    out = char(string(in));
    out = regexprep(out, '\\(it|rm|bf)\s*', '');
    out = strtrim(out);
end

function out = value_with_unit(val, unit, nfmt)
    sval = sigfig(val, nfmt);
    u = legend_plain_label(unit);
    if isempty(u)
        out = sval;
    else
        out = [sval ' ' u];
    end
end

function s = sigfig(val, n)
    if nargin < 2 || isempty(n), n = 4; end
    if isempty(val) || ~isfinite(val)
        s = num2str(val);
    else
        s = sprintf('%.*g', n, val);
    end
end


% Disable grid for cleaner uncertainty display
grid off;
