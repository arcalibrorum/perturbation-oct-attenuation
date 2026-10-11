function result = left_cumulative_integral(values, spacing)
%LEFT_CUMULATIVE_INTEGRAL Left-rectangle cumulative integral with result(1)=0.

values = validate_vector(values, 'values');
validateattributes(spacing, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
result = zeros(size(values));
result(2:end) = spacing .* cumsum(values(1:end-1));
end

function vector = validate_vector(values, name)
validateattributes(values, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'}, '', name);
if numel(values) < 2
    error('octattenuation:TooFewSamples', '%s must contain at least two samples.', name);
end
vector = double(values(:));
end
