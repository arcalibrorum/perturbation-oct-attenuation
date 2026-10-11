function result = reconstruct_piecewise(z, intensity, interval_edges, fixed_parameters, fit_windows, mu_bounds)
%RECONSTRUCT_PIECEWISE Apply independent local reconstruction by interval.
%   fixed_parameters is an intervals-by-2 array [mu0, amplitude]. If empty,
%   each interval is fitted. fit_windows is a cell array of local [start stop]
%   windows. A sample on an internal edge belongs to the deeper interval.

if nargin < 4
    fixed_parameters = [];
end
if nargin < 5
    fit_windows = [];
end
if nargin < 6 || isempty(mu_bounds)
    mu_bounds = [1e-5, 50];
end
validateattributes(z, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
validateattributes(intensity, {'numeric'}, {'vector', 'real', 'finite', 'nonnegative', 'nonempty'});
validateattributes(interval_edges, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
z = double(z(:));
intensity = double(intensity(:));
edges = double(interval_edges(:));
if numel(z) ~= numel(intensity) || numel(z) < 2 || any(diff(z) <= 0)
    error('octattenuation:InvalidInput', 'z and intensity must be equal-length vectors on an increasing grid.');
end
if numel(edges) < 2 || any(diff(edges) <= 0)
    error('octattenuation:InvalidEdges', 'interval_edges must be strictly increasing.');
end
interval_count = numel(edges) - 1;
if ~isempty(fixed_parameters)
    valid_numeric = isnumeric(fixed_parameters) && isequal(size(fixed_parameters), [interval_count, 2]);
    valid_struct = isstruct(fixed_parameters) && numel(fixed_parameters) == interval_count && ...
        all(isfield(fixed_parameters, {'mu0', 'amplitude'}));
    if ~valid_numeric && ~valid_struct
        error('octattenuation:InvalidParameters', ...
            'fixed_parameters must be an intervals-by-2 array or a fit structure per interval.');
    end
end
if ~isempty(fit_windows) && numel(fit_windows) ~= interval_count
    error('octattenuation:InvalidWindows', 'fit_windows must contain one entry per interval.');
end

mu0_map = nan(size(z));
mu1_map = nan(size(z));
fits = repmat(struct('mu0', NaN, 'amplitude', NaN, 'rmse', NaN, ...
    'r_squared', NaN, 'sample_count', 0), interval_count, 1);
for index = 1:interval_count
    start_index = find(z >= edges(index), 1, 'first');
    if index == interval_count
        stop_index = find(z <= edges(index + 1), 1, 'last');
    else
        stop_index = find(z < edges(index + 1), 1, 'last');
    end
    if isempty(start_index) || isempty(stop_index) || stop_index - start_index + 1 < 3
        error('octattenuation:TooFewSamples', 'Interval %d contains fewer than three samples.', index);
    end
    local_z = z(start_index:stop_index);
    local_i = intensity(start_index:stop_index);
    if isempty(fixed_parameters)
        if isempty(fit_windows)
            window = [];
        else
            window = fit_windows{index};
        end
        fit_result = octattenuation.fit_zeroth_order(local_z, local_i, mu_bounds, window);
    elseif isnumeric(fixed_parameters)
        fit_result = struct('mu0', fixed_parameters(index, 1), ...
            'amplitude', fixed_parameters(index, 2), 'rmse', NaN, ...
            'r_squared', NaN, 'sample_count', numel(local_z));
    else
        fit_result = fixed_parameters(index);
    end
    if ~isfield(fit_result, 'rmse'); fit_result.rmse = NaN; end
    if ~isfield(fit_result, 'r_squared'); fit_result.r_squared = NaN; end
    if ~isfield(fit_result, 'sample_count'); fit_result.sample_count = numel(local_z); end
    local = octattenuation.reconstruct_interval(local_z, local_i, fit_result);
    mu0_map(start_index:stop_index) = local.mu0;
    mu1_map(start_index:stop_index) = local.mu1;
    fits(index) = fit_result;
end
result = struct('mu0', mu0_map, 'mu1', mu1_map, 'total', mu0_map + mu1_map, 'fits', fits);
end
