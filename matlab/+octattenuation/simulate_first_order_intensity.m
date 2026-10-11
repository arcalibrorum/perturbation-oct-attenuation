function intensity = simulate_first_order_intensity(z, mu0, mu1, amplitude)
%SIMULATE_FIRST_ORDER_INTENSITY Generate the local first-order model signal.

if nargin < 4
    amplitude = 1;
end
validateattributes(z, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
validateattributes(mu1, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
validateattributes(mu0, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(amplitude, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
z = double(z(:));
mu1 = double(mu1(:));
if numel(z) ~= numel(mu1) || numel(z) < 2 || any(diff(z) <= 0)
    error('octattenuation:InvalidInput', 'z and mu1 must be equal-length vectors on an increasing grid.');
end
increments = diff(z);
spacing = median(increments);
if any(abs(increments - spacing) > max(1e-12, abs(spacing) * 2e-7))
    error('octattenuation:NonuniformGrid', 'A uniformly sampled depth grid is required.');
end
xi = z - z(1);
integral = octattenuation.left_cumulative_integral(mu1, spacing);
intensity = amplitude .* (mu0 + mu1 - 2 .* mu0 .* integral) .* exp(-2 .* mu0 .* xi);
if any(intensity < 0)
    error('octattenuation:NegativeIntensity', 'The correction produces negative first-order intensity.');
end
end
