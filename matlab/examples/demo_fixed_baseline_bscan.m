function demo_fixed_baseline_bscan(output_directory, make_plot)
%DEMO_FIXED_BASELINE_BSCAN Synthetic no-vessel-prior reconstruction workflow.

if nargin < 1 || isempty(output_directory)
    output_directory = fullfile(pwd, 'output', 'fixed_baseline_bscan_matlab');
end
if nargin < 2
    make_plot = true;
end
repository_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repository_root, 'matlab'));
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end
rng(7, 'twister');
z = (0:359).' .* 0.002;
x = linspace(-1, 1, 121);
upper = 0.08 + 0.018 .* (1 - x .^ 2);
middle = 0.34 + 0.025 .* (1 - x .^ 2);
lower = repmat(z(end), size(x));
boundaries = [upper; middle; lower];
parameters = [0.35, 2.0; 0.75, 1.2];
intensity = repmat(1e-5, numel(z), numel(x));
reference = nan(size(intensity));
vessel_mask = false(size(intensity));

for column = 1:numel(x)
    for layer = 1:2
        start_index = find(z >= boundaries(layer, column), 1, 'first');
        if layer == 2
            stop_index = find(z <= boundaries(layer + 1, column), 1, 'last');
        else
            stop_index = find(z < boundaries(layer + 1, column), 1, 'last');
        end
        local_z = z(start_index:stop_index);
        mu0 = parameters(layer, 1);
        amplitude = parameters(layer, 2);
        mu1 = zeros(size(local_z));
        if layer == 1
            vessel = exp(-(x(column) ./ 0.18) .^ 2) .* exp(-((local_z - 0.20) ./ 0.055) .^ 2);
            mu1 = 0.20 .* vessel;
            vessel_mask(start_index:stop_index, column) = vessel > 0.55;
        end
        intensity(start_index:stop_index, column) = ...
            octattenuation.simulate_first_order_intensity(local_z, mu0, mu1, amplitude);
        reference(start_index:stop_index, column) = mu0 + mu1;
    end
end
intensity = intensity .* max(0, 1 + 0.01 .* randn(size(intensity)));
control_columns = find(abs(x) > 0.65);
fitted_parameters = octattenuation.fit_control_baselines(z, intensity, boundaries, ...
    control_columns, [], [0.05, 2.0]);
reconstruction = octattenuation.reconstruct_bscan(z, intensity, boundaries, fitted_parameters);
surrounding_mask = false(size(vessel_mask));
depth_mask = any(vessel_mask, 2);
surrounding_columns = abs(x) > 0.30 & abs(x) < 0.55;
surrounding_mask(:, surrounding_columns) = repmat(depth_mask, 1, nnz(surrounding_columns));
[contrast, mae] = octattenuation.ablation_metrics(reconstruction.total, reference, ...
    vessel_mask, surrounding_mask);
save(fullfile(output_directory, 'fixed_baseline_case.mat'), 'z', 'x', 'intensity', ...
    'reference', 'vessel_mask', 'surrounding_mask', 'fitted_parameters', ...
    'reconstruction', 'contrast', 'mae');
fprintf('vessel contrast: %.6f mm^-1\n', contrast);
fprintf('vessel ROI MAE: %.6f mm^-1\n', mae);

if ~make_plot
    return;
end
figure_handle = figure('Color', 'w', 'Position', [100, 100, 1200, 300]);
tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
images = {intensity, reconstruction.mu0, abs(reconstruction.mu1), reconstruction.total};
titles = {'Linear intensity', '$\hat{\mu}_0$', '$|\hat{\mu}_1|$', '$\hat{\mu}_0+\hat{\mu}_1$'};
for panel = 1:4
    nexttile;
    imagesc(x, z, images{panel}); axis image; set(gca, 'YDir', 'reverse');
    title(titles{panel}, 'Interpreter', 'latex'); xlabel('Lateral coordinate'); colorbar;
end
exportgraphics(figure_handle, fullfile(output_directory, 'fixed_baseline_bscan_demo.png'), 'Resolution', 200);
end
