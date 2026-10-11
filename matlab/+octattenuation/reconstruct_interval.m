function result = reconstruct_interval(z, intensity, fit_result, mu_bounds, fit_window)
%RECONSTRUCT_INTERVAL Reconstruct one local depth interval.

if nargin < 3
    fit_result = [];
end
if nargin < 4 || isempty(mu_bounds)
    mu_bounds = [1e-5, 50];
end
if nargin < 5
    fit_window = [];
end
z = double(z(:));
intensity = double(intensity(:));
if isempty(fit_result)
    fit_result = octattenuation.fit_zeroth_order(z, intensity, mu_bounds, fit_window);
end
required = {'mu0', 'amplitude'};
if ~isstruct(fit_result) || ~all(isfield(fit_result, required))
    error('octattenuation:InvalidFit', 'fit_result must contain mu0 and amplitude.');
end
mu1 = octattenuation.first_order_correction(z, intensity, fit_result.mu0, fit_result.amplitude);
mu0 = repmat(double(fit_result.mu0), size(z));
result = struct('mu0', mu0, 'mu1', mu1, 'total', mu0 + mu1, 'fits', fit_result);
end
