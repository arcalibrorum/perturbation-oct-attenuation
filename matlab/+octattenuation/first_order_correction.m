function mu1 = first_order_correction(z, intensity, mu0, amplitude)
%FIRST_ORDER_CORRECTION Evaluate the manuscript's local first-order term.
%   z is in mm, intensity is linear, and mu0 is in mm^-1. The local
%   coordinate and integral are reset at the first sample.

[z, intensity] = validate_inputs(z, intensity);
validateattributes(mu0, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(amplitude, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
spacing = uniform_spacing(z);
xi = z - z(1);
cumulative_intensity = octattenuation.left_cumulative_integral(intensity, spacing);
bracket = intensity ./ amplitude + (2 .* mu0 ./ amplitude) .* cumulative_intensity - mu0;
mu1 = bracket .* exp(2 .* mu0 .* xi);
end

function spacing = uniform_spacing(z)
increments = diff(z);
spacing = median(increments);
tolerance = max(1e-12, abs(spacing) * 1e-10);
if any(abs(increments - spacing) > tolerance + 1e-7 * abs(spacing))
    error('octattenuation:NonuniformGrid', 'A uniformly sampled depth grid is required.');
end
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
