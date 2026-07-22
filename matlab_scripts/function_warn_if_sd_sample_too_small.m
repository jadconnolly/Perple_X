function factor = function_warn_if_sd_sample_too_small(label, n, target_rse, suggestion)
% Warn if the sample size is too small to estimate a standard deviation
% with relative standard error <= target_rse.
%
% Approximation used: RSE[s] ~ 1/sqrt(2*(n-1)).

if nargin < 3 || isempty(target_rse)
    target_rse = 0.20;
end
if nargin < 4
    suggestion = 'Increase the relevant run sizes proportionally';
end
assert(isscalar(target_rse) && isnumeric(target_rse) && isfinite(target_rse) && target_rse > 0 && target_rse < 1, ...
    'target_rse must be a finite scalar between 0 and 1')

n_required = ceil(1 + 1/(2*target_rse^2));

if isempty(n) || isnan(n) || n < 2
    factor = Inf;
    msg = sprintf(['Only %g usable point(s) were found in the %s. ', ...
        'At least %d are needed for an approximately %.0f%%-accurate standard deviation estimate. ', ...
        '%s.'], n, label, n_required, 100*target_rse, suggestion);
    warning(msg)
    try
        warndlg(msg, 'MC\_fit warning');
    catch
    end
    return
end

factor = n_required / n;
if n < n_required
    msg = sprintf(['The %s has only %d usable points. ', ...
        'To estimate a standard deviation to about %.0f%%, use about %d points ', ...
        '(increase by a factor of %.2f). %s by approximately that factor.'], ...
        label, n, 100*target_rse, n_required, factor, suggestion);
    warning(msg)
    try
        warndlg(msg, 'MC\_fit warning');
    catch
    end
end
end
