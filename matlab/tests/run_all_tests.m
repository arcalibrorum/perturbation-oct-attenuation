function run_all_tests
%RUN_ALL_TESTS Dependency-free numerical verification for the MATLAB package.

repository_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repository_root, 'matlab'));

z = (0:500).' .* 0.001;
expected_mu0 = 0.8;
expected_amplitude = 2.4;
intensity = expected_amplitude .* expected_mu0 .* exp(-2 .* expected_mu0 .* z);
fitted = octattenuation.fit_zeroth_order(z, intensity, [0.01, 3.0]);
assert(abs(fitted.mu0 - expected_mu0) < 1e-8, 'Zeroth-order attenuation fit failed.');
assert(abs(fitted.amplitude - expected_amplitude) < 1e-8, 'Zeroth-order amplitude fit failed.');

expected_mu1 = 0.05 .* sin(2 .* pi .* z ./ z(end));
first_order_intensity = octattenuation.simulate_first_order_intensity(z, expected_mu0, ...
    expected_mu1, expected_amplitude);
recovered_mu1 = octattenuation.first_order_correction(z, first_order_intensity, ...
    expected_mu0, expected_amplitude);
assert(max(abs(recovered_mu1 - expected_mu1)) < 8.1e-4, 'First-order inversion failed.');

simple = [4; 3; 2; 1];
expected_conventional = simple ./ (0.2 .* [10; 6; 3; 1]);
actual_conventional = octattenuation.conventional_depth_resolved(simple, 0.1);
assert(max(abs(actual_conventional - expected_conventional)) < 1e-12, 'Conventional estimator failed.');

control_z = (0:100).' .* 0.005;
control_boundaries = [0, 0.25, 0.5];
control_image = zeros(numel(control_z), 5);
upper_mask = control_z < 0.25;
for column = 1:size(control_image, 2)
    control_image(upper_mask, column) = 1.6 .* 0.4 .* exp(-2 .* 0.4 .* control_z(upper_mask));
    local_depth = control_z(~upper_mask) - control_z(find(~upper_mask, 1, 'first'));
    control_image(~upper_mask, column) = 1.1 .* 0.9 .* exp(-2 .* 0.9 .* local_depth);
end
control_image(:, 3) = 1.4 .* control_image(:, 3);
control_fits = octattenuation.fit_control_baselines(control_z, control_image, ...
    control_boundaries, [1, 2, 4, 5], [], [0.05, 2.0]);
assert(abs(control_fits(1).mu0 - 0.4) < 1e-7, 'First control baseline fit failed.');
assert(abs(control_fits(2).mu0 - 0.9) < 1e-7, 'Second control baseline fit failed.');
assert(control_fits(1).sample_count == 4, 'Control fit count is incorrect.');

shared_case = readtable(fullfile(repository_root, 'shared', 'test_data', 'cross_language_case.csv'));
shared_mu1 = octattenuation.first_order_correction(shared_case.z_mm, shared_case.intensity, 0.7, 1.8);
assert(max(abs(shared_mu1 - shared_case.mu1_expected)) < 2e-10, 'Cross-language case failed.');

fprintf('All MATLAB tests passed.\n');
end
