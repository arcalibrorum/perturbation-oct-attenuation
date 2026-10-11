function result = reconstruct_bscan(z, intensity, layer_boundaries, fixed_parameters)
%RECONSTRUCT_BSCAN Reconstruct a linear-intensity B-scan with fixed baselines.
%   intensity is depth-by-lateral. layer_boundaries may be a shared vector or
%   a (layers+1)-by-lateral array. No vessel mask is used by this routine.

validateattributes(z, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
validateattributes(intensity, {'numeric'}, {'2d', 'real', 'finite', 'nonnegative'});
z = double(z(:));
intensity = double(intensity);
if size(intensity, 1) ~= numel(z)
    error('octattenuation:InvalidShape', 'intensity must be depth-by-lateral.');
end
if isnumeric(fixed_parameters)
    validateattributes(fixed_parameters, {'numeric'}, {'2d', 'real', 'finite', 'positive'});
    if size(fixed_parameters, 2) ~= 2
        error('octattenuation:InvalidParameters', 'Numeric fixed_parameters must have two columns.');
    end
    layer_count = size(fixed_parameters, 1);
elseif isstruct(fixed_parameters) && all(isfield(fixed_parameters, {'mu0', 'amplitude'}))
    layer_count = numel(fixed_parameters);
else
    error('octattenuation:InvalidParameters', 'fixed_parameters must be numeric pairs or fit structures.');
end
if isvector(layer_boundaries)
    boundaries = repmat(double(layer_boundaries(:)), 1, size(intensity, 2));
else
    boundaries = double(layer_boundaries);
end
if ~isequal(size(boundaries), [layer_count + 1, size(intensity, 2)])
    error('octattenuation:InvalidBoundaries', 'layer_boundaries has an incompatible shape.');
end

mu0_map = nan(size(intensity));
mu1_map = nan(size(intensity));
for column = 1:size(intensity, 2)
    column_result = octattenuation.reconstruct_piecewise(z, intensity(:, column), ...
        boundaries(:, column), fixed_parameters);
    mu0_map(:, column) = column_result.mu0;
    mu1_map(:, column) = column_result.mu1;
end
result = struct('mu0', mu0_map, 'mu1', mu1_map, 'total', mu0_map + mu1_map, 'fits', []);
end
