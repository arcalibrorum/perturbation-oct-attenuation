function demo_finite_depth(output_directory, make_plot)
%DEMO_FINITE_DEPTH Public finite-depth A-line demonstration.

if nargin < 1 || isempty(output_directory)
    output_directory = fullfile(pwd, 'output', 'finite_depth_matlab');
end
if nargin < 2
    make_plot = true;
end
repository_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repository_root, 'matlab'));
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

sample_count = 1000;
spacing = 0.002;
z = (0:sample_count-1).' .* spacing;
step_mu = [0.6 .* ones(150, 1); 0.2 .* ones(10, 1); 0.6 .* ones(140, 1); ...
    0.4 .* ones(300, 1); 0.8 .* ones(100, 1); 1.0 .* ones(10, 1); ...
    0.8 .* ones(290, 1)];
smooth_mu = [ones(300, 1); 1.6 .* ones(300, 1); 1.2 .* ones(400, 1)] + ...
    0.05 .* sin(0.2 .* (1:1000).');
edges = [0, 0.6, 1.2, z(end)];
step_fit_windows = {[0, 0.298], [0, 0.598], [0, 0.198]};
smooth_fit_windows = {[0, 0.598], [0, 0.598], [0, 0.798]};

step_intensity = octattenuation.simulate_beer_lambert_intensity(z, step_mu, 1.67);
smooth_intensity = octattenuation.simulate_beer_lambert_intensity(z, smooth_mu, 1.0);
step_result = octattenuation.reconstruct_piecewise(z, step_intensity, edges, [], step_fit_windows);
smooth_result = octattenuation.reconstruct_piecewise(z, smooth_intensity, edges, [], smooth_fit_windows);
step_conventional = octattenuation.conventional_depth_resolved(step_intensity, spacing);
smooth_conventional = octattenuation.conventional_depth_resolved(smooth_intensity, spacing);

case_name = [repmat("step", sample_count, 1); repmat("smooth", sample_count, 1)];
table_result = table(case_name, [z; z], [step_intensity; smooth_intensity], ...
    [step_mu; smooth_mu], [step_result.mu0; smooth_result.mu0], ...
    [step_result.mu1; smooth_result.mu1], [step_result.total; smooth_result.total], ...
    [step_conventional; smooth_conventional], ...
    'VariableNames', {'case_name', 'z_mm', 'intensity', 'reference_mu', ...
    'mu0', 'mu1', 'total', 'conventional'});
writetable(table_result, fullfile(output_directory, 'finite_depth_results.csv'));

if ~make_plot
    return;
end
figure_handle = figure('Color', 'w', 'Position', [100, 100, 1000, 600]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; semilogy(z, step_intensity, 'Color', [0.25, 0.25, 0.25]);
title('Step-like linear OCT intensity'); xlabel('Depth (mm)'); ylabel('Intensity (a.u.)');
nexttile; semilogy(z, smooth_intensity, 'Color', [0.25, 0.25, 0.25]);
title('Smooth linear OCT intensity'); xlabel('Depth (mm)'); ylabel('Intensity (a.u.)');
nexttile; hold on;
plot(z, step_mu, 'k--'); plot(z, step_result.total); plot(z, step_conventional);
ylim([0, 3]); xlabel('Depth (mm)'); ylabel('Attenuation (mm^{-1})');
legend('Reference', 'Perturbation-based', 'Conventional', 'Location', 'best');
nexttile; hold on;
plot(z, smooth_mu, 'k--'); plot(z, smooth_result.total); plot(z, smooth_conventional);
ylim([0, 3]); xlabel('Depth (mm)'); ylabel('Attenuation (mm^{-1})');
exportgraphics(figure_handle, fullfile(output_directory, 'finite_depth_demo.png'), 'Resolution', 200);
end
