function result = conventional_depth_resolved(intensity, spacing, dimension)
%CONVENTIONAL_DEPTH_RESOLVED Conventional distal-tail attenuation estimator.

validateattributes(intensity, {'numeric'}, {'real', 'finite', 'nonnegative', 'nonempty'});
validateattributes(spacing, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
if nargin < 3 || isempty(dimension)
    dimension = find(size(intensity) ~= 1, 1, 'first');
    if isempty(dimension)
        dimension = 1;
    end
end
validateattributes(dimension, {'numeric'}, {'scalar', 'integer', 'positive'});
signal = double(intensity);
tail = flip(cumsum(flip(signal, dimension), dimension), dimension);
denominator = 2 .* spacing .* tail;
result = nan(size(signal));
valid = denominator > 0;
result(valid) = signal(valid) ./ denominator(valid);
end
