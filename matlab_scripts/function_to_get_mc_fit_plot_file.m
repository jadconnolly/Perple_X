function [x, y, z, symb, fit, xname, yname, zname, sname, uname, nrow, xvar, ...
    yvar, zvar, file, path, model, filled_symbols_for_nofit_auto, imprecision_filter_info] = ...
    function_to_get_mc_fit_plot_file( ...
    type, xvar, yvar, zvar, epsz, fitonly_prompt, manual_filter, filter, ...
    xname, yname, zname, sname, file, path, first, ...
    convert_to_MPa, convert_to_kbar, convert_to_GPa, convert_to_C, force_fit_filter, ...
    outlier_filter, outlier_filter_iterations, outlier_filter_convergence, ...
    outlier_filter_z, outlier_filter_delta_z, outlier_filter_delta_z_max)

% Read an mc_fit *.pts file and return the selected plotting variables.

filled_symbols_for_nofit_auto = false;
imprecision_filter_info = struct( ...
    'applied', false, ...
    'pre_count', NaN, ...
    'post_count', NaN, ...
    'removed_count', NaN, ...
    'epsz', NaN, ...
    'zbest', NaN);

if nargin < 20 || isempty(force_fit_filter)
    force_fit_filter = false;
end
if nargin < 21 || isempty(outlier_filter)
    outlier_filter = true;
end
if nargin < 22 || isempty(outlier_filter_iterations)
    outlier_filter_iterations = 10;
end
if nargin < 23 || isempty(outlier_filter_convergence)
    outlier_filter_convergence = 0.01;
end
if nargin < 24 || isempty(outlier_filter_z)
    outlier_filter_z = true;
end
if nargin < 25 || isempty(outlier_filter_delta_z)
    outlier_filter_delta_z = false;
end
if nargin < 26 || isempty(outlier_filter_delta_z_max)
    outlier_filter_delta_z_max = 3.91;
end

% Validate unit-conversion switches.
n_pressure_conversions = double(convert_to_MPa) + double(convert_to_kbar) + double(convert_to_GPa);
if n_pressure_conversions > 1
    error(['Only one pressure conversion may be true: ', ...
        'convert_to_MPa, convert_to_kbar, and convert_to_GPa are mutually exclusive.'])
end

bug = false;

if first && ~bug && (isempty(file) || isempty(path))
    [file, path, indx] = mc_fit_dialogs('select_plot_file', type);
    if indx == 0
        error(['You did not choose a ',type,' file, I quit!'])
    end
elseif first && bug
    file = 'NIL16_perturbed.pts';
    path = 'C:\jamie\Perple_X\ugh2\';
end

data_file = fullfile(path, file);
T = readtable(data_file, 'FileType','delimitedtext', 'PreserveVariableNames',true);

if isempty(T) || height(T) == 0 || width(T) == 0
    error(['File "', file, '" contains no data. ', ...
        'Verify that MC_fit completed successfully.'])
end

if ~isempty(regexp(file,'_central\.pts$','once','ignorecase'))
    text = 'the central';
    model = "Central";
elseif ~isempty(regexp(file,'_perturbed\.pts$','once','ignorecase'))
    text = 'the perturbed';
    model = "Perturbation";
else
    error(['No support for ordinary ',type,' files, I quit!'])
end

colnames = string(T.Properties.VariableNames);
A = table2array(T);
[nrow, cols] = size(A);

% Find end of parameter columns (999 marker lies after the last parameter).
cols = find(A(1,1:cols)==999, 1, 'first') - 1;

bayes = 0;
if any(A(:,1)==5) || any(A(:,1)==6)
    bayes = 3;
end

if bayes == 0
    zname = "ln(Misfit)";
else
    zname = "ln(Bayes)";
end

% Unit conversions.
if convert_to_MPa
    A(:,3) = A(:,3)/10;        % P: bar -> MPa
    pressure_unit = "MPa";
elseif convert_to_GPa
    A(:,3) = A(:,3)/10000;     % P: bar -> GPa
    pressure_unit = "GPa";
elseif convert_to_kbar
    A(:,3) = A(:,3)/1000;      % P: bar -> kbar
    pressure_unit = "kbar";
else
    pressure_unit = "bar";
end

if convert_to_C
    A(:,4) = A(:,4) - 273.15;  % T: K -> °C
    temperature_unit = '\circC';
else
    temperature_unit = "K";
end

A(:,cols) = log(A(:,cols));

% Friendly variable names for selection dialog.
% Variables available for x/y plotting are columns 3:(cols-1).
% The final column (OBJF) is reserved for z and is not selectable.
xyvars = 3:(cols-1);
dnames = cell(1, numel(xyvars));
for j = 1:numel(xyvars)
    i = xyvars(j);
    if i == 3
        dnames{j} = ['P, ' char(pressure_unit)];
    elseif i == 4
        if convert_to_C
            dnames{j} = 'T, C';
        else
            dnames{j} = ['T, ' char(temperature_unit)];
        end
    else
        dnames{j} = char(string(colnames(i)));
    end
end

if first
    if ~bug
        [xpick, ok] = mc_fit_dialogs('select_axis_variable', 'x', dnames);
        if ok == 0, error('You did not choose a variable, I quit!'), end

        ymask = true(1, numel(dnames));
        ymask(xpick) = false;
        ynames = dnames(ymask);
        ymap = find(ymask);

        [ypick_local, ok] = mc_fit_dialogs('select_axis_variable', 'y', ynames);
        if ok == 0, error('You did not choose a variable, I quit!'), end
        ypick = ymap(ypick_local);
    else
        xpick = min(2, numel(xyvars));
        ypick = 1;
    end

    xvar = xyvars(xpick);
    yvar = xyvars(ypick);
    zvar = cols;
end

[xname, sname(1), uname(1)] = local_make_axis_label(xvar, colnames(xvar), pressure_unit, temperature_unit);
[yname, sname(2), uname(2)] = local_make_axis_label(yvar, colnames(yvar), pressure_unit, temperature_unit);
sname(3) = zname;
uname(3) = "";

symb = A(:,1);
fit  = A(:,2);
x    = A(:,xvar);
y    = A(:,yvar);
z    = A(:,zvar);

% Identify the best central model and suppress any ordinary-trial replicate
% immediately after reading.  The best central model may appear twice in
% central *.pts files: once as the explicit best-model marker (symb==1) and
% once as an ordinary central trial (symb==3, or symb==6 for Bayes output).
% The ordinary-trial replicate is an output artifact, not an independent
% model, and must not enter fit-all-data counts, imprecision filtering,
% retained-fraction calculations, covariance estimates, or warning logic.
% If the explicit symb==1 marker is absent, promote the best finite central
% trial to symb==1 before suppressing any remaining replicate.
if model == "Central"
    cent_sym = 3 + 3*double(zname == "ln(Bayes)");
    [symb, x, y, z] = local_prepare_best_central_model(symb, x, y, z, cent_sym, zname);
end

% Optional robust outlier filtering for perturbation models. This is done
% immediately after the objective function has been transformed to
% ln(Misfit)/ln(Bayes), so all downstream perturbation-derived imprecision
% estimates, and the central-model filtering based on those estimates, use the
% filtered perturbation cloud. The central-model data are not outlier-filtered.

if outlier_filter && model == "Perturbation"
    pert_sym = 2 + 3*double(zname == "ln(Bayes)");
    idx_out = [];

    if outlier_filter_z && outlier_filter_iterations > 0
        idx_z = local_robust_one_sided_z_outliers( ...
            z, symb, pert_sym, zname, outlier_filter_iterations, outlier_filter_convergence);
        idx_out = union(idx_out, idx_z);
    end

    if outlier_filter_delta_z
        z_for_delta = z;
        z_for_delta(idx_out) = NaN;   % report only additional delta-z eliminations
        idx_dz = local_delta_z_outliers( ...
            z_for_delta, symb, pert_sym, zname, outlier_filter_delta_z_max);
        idx_out = union(idx_out, idx_dz);
    end

    if ~isempty(idx_out)
        % Safety guard: outlier filtering is only allowed to suppress
        % perturbation trial rows. The perturbation *.pts file also carries
        % the best central model marker (symb==1), which must remain finite
        % so downstream plotting can identify and draw the best central model.
        idx_out = idx_out(symb(idx_out) == pert_sym);
        if ~isempty(idx_out)
            x(idx_out) = NaN;
            y(idx_out) = NaN;
            z(idx_out) = NaN;
        end
    end
end

if model == "Central"
    % Include the best central model marker (symb==1) so it is removed too
    % when the user chooses to retain only fitting central models.
    istrial = ismember(symb, [1 3 6]);
else
    istrial = ismember(symb, [2 5]);
end
sall  = sum(istrial);
sins  = sum(istrial & fit==1);
souts = sum(istrial & fit==0);

% Automatic marker convention: use filled symbols for non-fitting models only
% when no trial model fits all data. If all models fit, or there is a mixture
% of fitting and non-fitting models, use the opposite convention.
filled_symbols_for_nofit_auto = (sall > 0 && sall == souts);

if fitonly_prompt || force_fit_filter
    if sall ~= sins && sall ~= souts
        if force_fit_filter
            do_filter = true;
        else
            if model == "Central"
                model_text = 'central';
            else
                model_text = 'perturbation';
            end
            do_filter = mc_fit_dialogs('ask_fit_filter_single', sins, sall, model_text);
        end

        if do_filter
            n_considered = sum(istrial & isfinite(x) & isfinite(y) & isfinite(z));
            idx = find(istrial & fit==0);
            x(idx) = NaN; y(idx) = NaN; z(idx) = NaN;
            n_retained = sum(istrial & isfinite(x) & isfinite(y) & isfinite(z));
            fprintf('MC_fit fit-all-data filter (%s): considered %d models, retained %d, eliminated %d.\n', ...
                char(model), n_considered, n_retained, n_considered - n_retained);
        end
    end
end

if (filter && model == "Central") || manual_filter
    if manual_filter
        do_filter = mc_fit_dialogs('ask_manual_misfit_filter', text);
        if do_filter
            [lb_in, ub_in] = mc_fit_dialogs('ask_manual_misfit_limits', text, zname, ...
                min(z(:),[],'omitnan'), max(z(:),[],'omitnan'));
            if isempty(lb_in) || isempty(ub_in)
                return
            end
            lb = 0.999*lb_in;
            ub = 1.001*ub_in;
            finite_before = isfinite(x) & isfinite(y) & isfinite(z);
            n_considered = sum(finite_before);
            idx = find((z(:) < lb | z(:) > ub) & finite_before(:));
            x(idx) = NaN; y(idx) = NaN; z(idx) = NaN;
            n_retained = sum(isfinite(x) & isfinite(y) & isfinite(z));
            fprintf('MC_fit manual %s-value filter (%s): considered %d models, retained %d, eliminated %d.\n', ...
                char(zname), char(model), n_considered, n_retained, n_considered - n_retained);
        end
    else
        % Automatic imprecision filter for central models.  Count the finite
        % central-model sample immediately before and after this filter so
        % the plotting routine can distinguish a genuinely undersampled
        % central run from a well sampled run whose aleatoric scatter is
        % negligible relative to the perturbation-derived data uncertainty.
        central_trial = ismember(symb, [1 3 6]);
        finite_central_before = central_trial & isfinite(x) & isfinite(y) & isfinite(z);
        n_central_before = sum(finite_central_before);

        if zname == "ln(Bayes)"
            zbest = max(z(:),[],'omitnan');
        else
            zbest = min(z(:),[],'omitnan');
        end

        if isempty(epsz) || isnan(epsz)
            error('Filtering the central file requires a finite perturbation-derived z tolerance.')
        end

        lb = zbest - epsz;
        ub = zbest + epsz;
        idx = find(z(:) < lb | z(:) > ub);
        x(idx) = NaN; y(idx) = NaN; z(idx) = NaN;

        finite_central_after = central_trial & isfinite(x) & isfinite(y) & isfinite(z);
        n_central_after = sum(finite_central_after);

        imprecision_filter_info.applied = true;
        imprecision_filter_info.pre_count = n_central_before;
        imprecision_filter_info.post_count = n_central_after;
        imprecision_filter_info.removed_count = n_central_before - n_central_after;
        imprecision_filter_info.epsz = epsz;
        imprecision_filter_info.zbest = zbest;

        fprintf('MC_fit central imprecision filter: considered %d models, retained %d, eliminated %d.\n', ...
            n_central_before, n_central_after, n_central_before - n_central_after);
    end
end

end


function [symb, x, y, z] = local_prepare_best_central_model(symb, x, y, z, cent_sym, zname)
    idx_best = find(symb == 1 & isfinite(x) & isfinite(y) & isfinite(z));
    idx_cent = find(symb == cent_sym & isfinite(x) & isfinite(y) & isfinite(z));

    if isempty(idx_best) && ~isempty(idx_cent)
        if zname == "ln(Bayes)"
            [~, k] = max(z(idx_cent));
        else
            [~, k] = min(z(idx_cent));
        end
        idx_best = idx_cent(k);
        symb(idx_best) = 1;
        idx_cent(idx_cent == idx_best) = [];
    elseif ~isempty(idx_best)
        idx_best = idx_best(1);
    else
        return
    end

    if isempty(idx_cent)
        return
    end

    idx_rep = local_find_best_central_replicates(idx_cent, x, y, z, idx_best);
    if ~isempty(idx_rep)
        x(idx_rep) = NaN;
        y(idx_rep) = NaN;
        z(idx_rep) = NaN;
    end
end

function idx_rep = local_find_best_central_replicates(idx_cent, x, y, z, idx_best)
    idx_rep = [];
    if isempty(idx_cent) || isempty(idx_best)
        return
    end

    xb = x(idx_best);
    yb = y(idx_best);
    zb = z(idx_best);

    sx = max(1, max(abs([x(idx_cent); xb]), [], 'omitnan'));
    sy = max(1, max(abs([y(idx_cent); yb]), [], 'omitnan'));
    sz = max(1, max(abs([z(idx_cent); zb]), [], 'omitnan'));
    tol = 100*eps;

    is_rep = abs(x(idx_cent) - xb) <= tol*sx & ...
             abs(y(idx_cent) - yb) <= tol*sy & ...
             abs(z(idx_cent) - zb) <= tol*sz;

    idx_rep = idx_cent(is_rep);
end

function idx_out = local_robust_one_sided_z_outliers(z, symb, pert_sym, zname, niter, convergence_tol)
    idx_out = [];
    if nargin < 6 || isempty(convergence_tol)
        convergence_tol = 0.01;
    end

    idx = find(symb == pert_sym & isfinite(z));
    if numel(idx) < 3
        fprintf('MC_fit perturbation robust one-sided z filter: considered %d models, retained %d, eliminated 0 (fewer than 3 candidates).\n', ...
            numel(idx), numel(idx));
        return
    end

    % Robust one-sided z filter on the objective-function coordinate.
    % The old implementation used a mean/variance squared z-score
    % ("1-D Mahalanobis") and rejected both high and low tails. That was
    % vulnerable to bad-model clusters inflating the variance, and it could
    % also reject unusually good models. This version centers on the median,
    % scales by 1.4826*MAD for normal-equivalent sigma units, and rejects
    % only the bad tail:
    %   ln(Misfit): reject values above the median by more than the cutoff
    %   ln(Bayes):  reject values below the median by more than the cutoff
    % The best retained perturbation model is explicitly protected.
    robust_z_cutoff = sqrt(6.63489660102121);  % 99% normal-equivalent cutoff, avoids Statistics Toolbox dependency.

    keep = true(size(idx));

    for iter = 1:niter
        j = idx(keep);
        if numel(j) < 3
            fprintf('MC_fit perturbation robust one-sided z filter iteration %d: considered %d models, retained %d, eliminated 0 (fewer than 3 retained).\n', ...
                iter, numel(j), numel(j));
            break
        end

        zj = z(j);
        ctr = median(zj);
        mad_raw = median(abs(zj - ctr));
        scale = 1.4826 * mad_raw;

        if ~isfinite(scale) || scale <= 0
            fprintf('MC_fit perturbation robust one-sided z filter iteration %d: considered %d models, retained %d, eliminated 0 (zero/non-finite MAD scale).\n', ...
                iter, numel(j), numel(j));
            break
        end

        if zname == "ln(Bayes)"
            % For Bayes, larger is better, so reject the low-z tail only.
            robust_z = (ctr - zj) ./ scale;
            [~, kbest] = max(zj);
        else
            % For misfit, smaller is better, so reject the high-z tail only.
            robust_z = (zj - ctr) ./ scale;
            [~, kbest] = min(zj);
        end

        reject_local = robust_z > robust_z_cutoff;
        reject_local(kbest) = false;  % never reject the best retained perturbation model

        n_reject = sum(reject_local);
        kept_positions = find(keep);
        if n_reject > 0
            keep(kept_positions(reject_local)) = false;
        end

        fprintf('MC_fit perturbation robust one-sided z filter iteration %d: considered %d models, retained %d, eliminated %d.\n', ...
            iter, numel(j), numel(j) - n_reject, n_reject);

        if n_reject == 0
            break
        end

        if n_reject/numel(j) < convergence_tol
            break
        end
    end

    idx_out = idx(~keep);
end

function idx_out = local_delta_z_outliers(z, symb, pert_sym, zname, dzmax)
    idx_out = [];
    if isempty(dzmax) || ~isfinite(dzmax) || dzmax < 0
        fprintf('MC_fit perturbation delta-z filter: considered 0 models, retained 0, eliminated 0 (invalid cutoff).\n');
        return
    end

    idx = find(symb == pert_sym & isfinite(z));
    if isempty(idx)
        fprintf('MC_fit perturbation delta-z filter: considered 0 models, retained 0, eliminated 0 (no finite candidates).\n');
        return
    end

    if zname == "ln(Bayes)"
        % For Bayes, larger is better; reject models more than dzmax below best.
        zbest = max(z(idx), [], 'omitnan');
        dz = zbest - z(idx);
    else
        % For Misfit, smaller is better; reject models more than dzmax above best.
        zbest = min(z(idx), [], 'omitnan');
        dz = z(idx) - zbest;
    end

    idx_out = idx(dz > dzmax);
    fprintf('MC_fit perturbation delta-z filter: considered %d models, retained %d, eliminated %d (threshold %.4g in %s).\n', ...
        numel(idx), numel(idx) - numel(idx_out), numel(idx_out), dzmax, char(zname));
end

function [axis_name, sname, uname] = local_make_axis_label(var_index, raw_name, pressure_unit, temperature_unit)
    raw_name = char(string(raw_name));
    uname = "";
    if var_index == 4
        axis_name = "\itT\rm, " + temperature_unit;
        sname = "T";
        uname = " " + temperature_unit;
    elseif var_index == 3
        axis_name = "\itP\rm, " + pressure_unit;
        sname = "P";
        uname = " " + pressure_unit;
    else
        sname = string(raw_name);
        axis_name = local_format_symbol(string(raw_name));
    end
end

function out = local_format_symbol(name)
    name = string(name);
    if contains(name, "_")
        parts = split(name, "_");
        head = parts(1);
        tail = join(parts(2:end), "_");
        out = "\it" + head + "_{\rm" + tail + "}";
    else
        out = "\it" + name;
    end
end
