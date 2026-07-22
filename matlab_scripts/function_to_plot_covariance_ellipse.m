function [hcov, htot, avg, covariance, hpcov] = ...
    function_to_plot_covariance_ellipse(data, z, jplot, lw, pcovar, tot, sigma_level, center_xy)

% Plot covariance ellipse(s).
% Backward-compatible defaults are provided so older callers that only pass
% (data, z, jplot) or omit later options still work.

if nargin < 4 || isempty(lw)
    lw = 1.0;
end
if nargin < 5
    pcovar = [];
end
if nargin < 6 || isempty(tot)
    tot = false;
end
if nargin < 7 || isempty(sigma_level)
    sigma_level = 1;
end
if nargin < 8
    center_xy = [];
end

assert(isscalar(sigma_level) && isnumeric(sigma_level) && isfinite(sigma_level) && sigma_level > 0, ...
    'sigma_level must be a finite positive scalar')
assert(isempty(center_xy) || (isnumeric(center_xy) && numel(center_xy) == 2 && all(isfinite(center_xy(:)))), ...
    'center_xy must be empty or a finite two-element numeric vector')

% - data covariance in requested color
% - pcovar in red if tot == true
% - total covariance (covariance + pcovar) in black if tot == true

% --- base covariance ---
covariance = cov(data,"omitrows");

% Mean of the data is retained as an output for backward compatibility, but
% the plotted ellipse may be centered explicitly (e.g., on the best central
% model) rather than on the sample mean.
avg = mean(data,"omitnan");
if isempty(center_xy)
    X0 = avg(1);
    Y0 = avg(2);
else
    center_xy = center_xy(:).';
    X0 = center_xy(1);
    Y0 = center_xy(2);
end

% color selection
if jplot == 3
    lc = '-b';
elseif jplot == 2
    lc = '-r';
else
    lc = '-k';
end

% helper to compute ellipse from covariance matrix
    function r_ellipse = make_ellipse(C)
        [V,D] = eig(C);

        % largest eigenvalue/vector
        [largest_eigenval, idx] = max(diag(D));
        largest_eigenvec = V(:,idx);

        % smallest eigenvalue
        smallest_eigenval = min(diag(D));

        % angle
        phi = atan2(largest_eigenvec(2), largest_eigenvec(1));
        if phi < 0
            phi = phi + 2*pi;
        end

        % ellipse parameters

        p = erf(sigma_level / sqrt(2));   % map sigma → probability
        c = -2 * log(1 - p);              % 2D chi-square threshold
        chisquare_val = sqrt(c);
        theta = linspace(0,2*pi);

        a = chisquare_val * sqrt(largest_eigenval);
        b = chisquare_val * sqrt(smallest_eigenval);

        xr = a * cos(theta);
        yr = b * sin(theta);

        % rotation
        R = [cos(phi) sin(phi); -sin(phi) cos(phi)];
        r_ellipse = [xr; yr]' * R;
    end

% --- plot primary covariance ---
rell = make_ellipse(covariance);
zel = z * ones(size(rell,1),1);

hcov = plot3(rell(:,1)+X0, rell(:,2)+Y0, zel, lc, ...
    'LineWidth', lw);
hold on;

% --- initialize optional outputs ---
hpcov = [];
htot  = [];

% --- plot pcovar and total covariance if requested ---
if tot && ~isempty(pcovar)
    % plot previous covariance in red
    rell_p = make_ellipse(pcovar);
    hpcov = plot3(rell_p(:,1)+X0, rell_p(:,2)+Y0, zel, '-r', ...
        'LineWidth', lw);

    % plot total covariance in black
    Ctot = covariance + pcovar;   % independence assumption
    rell_tot = make_ellipse(Ctot);
    htot = plot3(rell_tot(:,1)+X0, rell_tot(:,2)+Y0, zel, '-k', ...
        'LineWidth', lw);
end

end
