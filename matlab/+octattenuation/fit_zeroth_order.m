function fit_result = fit_zeroth_order(z, intensity, mu_bounds, fit_window)
%FIT_ZEROTH_ORDER Fit I=A*mu0*exp(-2*mu0*xi) by least squares.
%   The amplitude is eliminated analytically for every candidate mu0. A
%   dependency-free golden-section search solves the remaining bounded
%   scalar minimization.

if nargin < 3 || isempty(mu_bounds)
    mu_bounds = [1e-5, 50];
end
if nargin < 4
    fit_window = [];
end
[z, intensity] = validate_inputs(z, intensity);
validateattributes(mu_bounds, {'numeric'}, {'vector', 'numel', 2, 'real', 'finite', 'positive'});
mu_bounds = double(mu_bounds(:)).';
if mu_bounds(1) >= mu_bounds(2)
    error('octattenuation:InvalidBounds', 'mu_bounds must satisfy lower < upper.');
end

xi = z - z(1);
mask = true(size(z));
if ~isempty(fit_window)
    validateattributes(fit_window, {'numeric'}, {'vector', 'numel', 2, 'real', 'finite'});
    if fit_window(1) >= fit_window(2)
        error('octattenuation:InvalidWindow', 'fit_window must satisfy start < stop.');
    end
    mask = xi >= fit_window(1) & xi <= fit_window(2);
end
if nnz(mask) < 3
    error('octattenuation:TooFewFitSamples', 'The fitting window must contain at least three samples.');
end
x_fit = xi(mask);
y_fit = intensity(mask);

objective = @(mu0) local_sse(mu0, x_fit, y_fit);
mu0 = golden_section(objective, mu_bounds(1), mu_bounds(2), 1e-11, 1000);
amplitude = local_amplitude(mu0, x_fit, y_fit);
prediction = amplitude .* mu0 .* exp(-2 .* mu0 .* x_fit);
residual = y_fit - prediction;
sse = sum(residual .^ 2);
rmse = sqrt(sse ./ numel(y_fit));
sst = sum((y_fit - mean(y_fit)) .^ 2);
if sst > 0
    r_squared = 1 - sse ./ sst;
else
    r_squared = NaN;
end
fit_result = struct('mu0', mu0, 'amplitude', amplitude, 'rmse', rmse, ...
    'r_squared', r_squared, 'sample_count', numel(y_fit));
end

function value = local_sse(mu0, x, y)
amplitude = local_amplitude(mu0, x, y);
residual = y - amplitude .* mu0 .* exp(-2 .* mu0 .* x);
value = sum(residual .^ 2);
end

function amplitude = local_amplitude(mu0, x, y)
basis = mu0 .* exp(-2 .* mu0 .* x);
denominator = sum(basis .^ 2);
if denominator <= realmin
    amplitude = 0;
else
    amplitude = max(0, sum(basis .* y) ./ denominator);
end
end

function optimum = golden_section(objective, lower, upper, tolerance, max_iterations)
ratio = (sqrt(5) - 1) / 2;
left = lower;
right = upper;
x1 = right - ratio * (right - left);
x2 = left + ratio * (right - left);
f1 = objective(x1);
f2 = objective(x2);
for iteration = 1:max_iterations %#ok<NASGU>
    if right - left <= tolerance * max([1, abs(x1), abs(x2)])
        break;
    end
    if f1 <= f2
        right = x2;
        x2 = x1;
        f2 = f1;
        x1 = right - ratio * (right - left);
        f1 = objective(x1);
    else
        left = x1;
        x1 = x2;
        f1 = f2;
        x2 = left + ratio * (right - left);
        f2 = objective(x2);
    end
end
optimum = (left + right) / 2;
end

function [z, intensity] = validate_inputs(z, intensity)
validateattributes(z, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'}, '', 'z');
validateattributes(intensity, {'numeric'}, {'vector', 'real', 'finite', 'nonnegative', 'nonempty'}, '', 'intensity');
z = double(z(:));
intensity = double(intensity(:));
if numel(z) < 2 || numel(z) ~= numel(intensity)
    error('octattenuation:InvalidShape', 'z and intensity must have the same length of at least two.');
end
if any(diff(z) <= 0)
    error('octattenuation:InvalidDepth', 'z must be strictly increasing.');
end
end
