function [contrast, mae] = ablation_metrics(image, reference, vessel_mask, background_mask)
%ABLATION_METRICS Vessel contrast and vessel-ROI mean absolute error.

validateattributes(image, {'numeric'}, {'real'});
validateattributes(reference, {'numeric'}, {'real', 'size', size(image)});
validateattributes(vessel_mask, {'logical', 'numeric'}, {'size', size(image)});
validateattributes(background_mask, {'logical', 'numeric'}, {'size', size(image)});
vessel_mask = logical(vessel_mask);
background_mask = logical(background_mask);
if ~any(vessel_mask, 'all') || ~any(background_mask, 'all')
    error('octattenuation:EmptyMask', 'Both masks must contain samples.');
end
contrast = median(image(vessel_mask), 'omitnan') - median(image(background_mask), 'omitnan');
mae = mean(abs(image(vessel_mask) - reference(vessel_mask)), 'omitnan');
end
