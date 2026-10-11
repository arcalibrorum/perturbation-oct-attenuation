function aggregate = fit_control_baselines(z, intensity, layer_boundaries, control_columns, fit_windows, mu_bounds)
%FIT_CONTROL_BASELINES Median layer fits from remote control A-lines.
%   This function does not accept or use a vessel mask. sample_count records
%   the number of accepted control-column fits for each layer.

if nargin < 5
    fit_windows = [];
end
if nargin < 6 || isempty(mu_bounds)
    mu_bounds = [1e-5, 50];
end
validateattributes(z, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
validateattributes(intensity, {'numeric'}, {'2d', 'real', 'finite', 'nonnegative'});
z = double(z(:));
intensity = double(intensity);
if size(intensity, 1) ~= numel(z)
    error('octattenuation:InvalidShape', 'intensity must be depth-by-lateral.');
end
if isvector(layer_boundaries)
    boundaries = repmat(double(layer_boundaries(:)), 1, size(intensity, 2));
else
    boundaries = double(layer_boundaries);
end
layer_count = size(boundaries, 1) - 1;
if ~isequal(size(boundaries), [layer_count + 1, size(intensity, 2)])
    error('octattenuation:InvalidBoundaries', 'layer_boundaries has an incompatible shape.');
end
if islogical(control_columns)
    if numel(control_columns) ~= size(intensity, 2)
        error('octattenuation:InvalidControls', 'The control mask must match the lateral dimension.');
    end
    columns = find(control_columns);
else
    columns = double(control_columns(:));
end
if isempty(columns) || any(columns < 1) || any(columns > size(intensity, 2)) || any(mod(columns, 1) ~= 0)
    error('octattenuation:InvalidControls', 'control_columns must contain valid MATLAB column indices.');
end
if ~isempty(fit_windows) && numel(fit_windows) ~= layer_count
    error('octattenuation:InvalidWindows', 'fit_windows must contain one entry per layer.');
end

aggregate = repmat(struct('mu0', NaN, 'amplitude', NaN, 'rmse', NaN, ...
    'r_squared', NaN, 'sample_count', 0), layer_count, 1);
for layer = 1:layer_count
    layer_fits = repmat(aggregate(1), numel(columns), 1);
    accepted = false(numel(columns), 1);
    for item = 1:numel(columns)
        column = columns(item);
        start_index = find(z >= boundaries(layer, column), 1, 'first');
        if layer == layer_count
            stop_index = find(z <= boundaries(layer + 1, column), 1, 'last');
        else
            stop_index = find(z < boundaries(layer + 1, column), 1, 'last');
        end
        if isempty(start_index) || isempty(stop_index) || stop_index - start_index + 1 < 3
            continue;
        end
        if isempty(fit_windows)
            window = [];
        else
            window = fit_windows{layer};
        end
        try
            fitted = octattenuation.fit_zeroth_order(z(start_index:stop_index), ...
                intensity(start_index:stop_index, column), mu_bounds, window);
        catch exception
            if startsWith(exception.identifier, 'octattenuation:')
                continue;
            end
            rethrow(exception);
        end
        if fitted.mu0 > 0 && fitted.amplitude > 0 && isfinite(fitted.mu0 + fitted.amplitude)
            layer_fits(item) = fitted;
            accepted(item) = true;
        end
    end
    if ~any(accepted)
        error('octattenuation:NoValidFit', 'No valid control-column fit was obtained for layer %d.', layer);
    end
    valid = layer_fits(accepted);
    aggregate(layer) = struct( ...
        'mu0', median([valid.mu0]), ...
        'amplitude', median([valid.amplitude]), ...
        'rmse', median([valid.rmse]), ...
        'r_squared', median([valid.r_squared]), ...
        'sample_count', numel(valid));
end
end
