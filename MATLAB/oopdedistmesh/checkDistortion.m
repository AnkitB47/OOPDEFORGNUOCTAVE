% distortion criteria function
function meets_criteria = checkDistortion(coordinates, elements)
    % Compute aspect ratio for each element
    aspect_ratios = computeAspectRatios(coordinates, elements);

    % Define a threshold for aspect ratio
    threshold_aspect_ratio = 2; % Adjust as needed

    % Check if all elements meet the aspect ratio criterion
    meets_criteria = all(aspect_ratios <= threshold_aspect_ratio);
end
