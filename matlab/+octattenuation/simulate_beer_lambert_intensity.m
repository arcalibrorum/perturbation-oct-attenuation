function intensity = simulate_beer_lambert_intensity(z, attenuation, amplitude)
%SIMULATE_BEER_LAMBERT_INTENSITY Ideal linear OCT intensity generator.

if nargin < 3
    amplitude = 1;
end
validateattributes(z, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
validateattributes(attenuation, {'numeric'}, {'vector', 'real', 'finite', 'nonnegative', 'nonempty'});
validateattributes(amplitude, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
z = double(z(:));
attenuation = double(attenuation(:));
if numel(z) ~= numel(attenuation) || numel(z) < 2 || any(diff(z) <= 0)
    error('octattenuation:InvalidInput', 'z and attenuation must be equal-length vectors on an increasing grid.');
end
spacing = uniform_spacing(z);
optical_depth = octattenuation.left_cumulative_integral(attenuation, spacing);
intensity = amplitude .* attenuation .* exp(-2 .* optical_depth);
end

function spacing = uniform_spacing(z)
increments = diff(z);
spacing = median(increments);
if any(abs(increments - spacing) > max(1e-12, abs(spacing) * 2e-7))
    error('octattenuation:NonuniformGrid', 'A uniformly sampled depth grid is required.');
end
end
