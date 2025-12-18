function overlay = Overlay_DGV(x, y, sizeArr, position, value, background, SNR, ifSNR, scale)
    % Overlay heatmap onto a background image
    if nargin < 7, SNR = []; end
    if nargin < 8, ifSNR = false; end

    heatmap = zeros(y, x, 'single');
    count_map = zeros(y, x, 'single');

    % Update heatmap for each point
    for i = 1:size(position, 1)
        x_center = round(position(i,1));
        y_center = round(position(i,2));
        score = value(i);
        if ifSNR
            snr = SNR(i);
        else
            snr = 1;
        end
        r = sizeArr(i);

        [heatmap, count_map] = update_heatmap(heatmap, x_center, y_center, score, snr, count_map, ifSNR, 1, r);
    end

    % Normalize heatmap
    mask = (count_map ~= 0);
    heatmap(mask) = heatmap(mask) ./ count_map(mask);

    % Clip and normalize
    vmin = scale(1); vmax = scale(2);
    heatmap_clipped = min(max(heatmap, vmin), vmax);
    heatmap_normalized = (heatmap_clipped - vmin) ./ (vmax - vmin);

    % Apply colormap (rainbow-like colormap adjusted)
    cmap = jet(256); % MATLABes parula or jet colormap
    idx = max(1, min(256, round(heatmap_normalized * 255) + 1));
    heatmap_colored = ind2rgb(idx, cmap);

    % Convert to uint8
    heatmap_colored = uint8(heatmap_colored * 255);

    % Transparency
    alpha = zeros(size(heatmap));
    alpha(heatmap > 0) = 0.4;

    overlay = background;
    for c = 1:3
        overlay(:,:,c) = uint8(double(background(:,:,c)) .* (1 - alpha) + double(heatmap_colored(:,:,c)) .* alpha);
    end
end

function val = gaussian_distance(d, sigma)
    % Gaussian function for distance weights
    val = exp(-(d.^2) / (2 * sigma^2));
end

function [heatmap, count_map] = update_heatmap(heatmap, x_center, y_center, score, max_snr, count_map, ifSNR, P, radius)
    if nargin < 7, ifSNR = true; end
    if nargin < 8, P = 1; end
    if nargin < 9, radius = 15; end

    [h, w] = size(heatmap);

    score(isnan(score))=0;

    % Define region of interest
    x_range = max(1, x_center-radius) : min(w, x_center+radius);
    y_range = max(1, y_center-radius) : min(h, y_center+radius);

    [x_grid, y_grid] = meshgrid(x_range, y_range);

    % Distances
    distances = sqrt((x_grid - x_center).^2 + (y_grid - y_center).^2);
    valid_mask = distances < radius;
    distances = distances(valid_mask);

    % Gaussian distance weights
    sigma = radius / 5;
    dist_weights = gaussian_distance(distances, sigma);

    % Apply SNR if needed
    if ifSNR
        snr_weights = max_snr * ones(size(dist_weights));
        weights = dist_weights .* (snr_weights.^P);
    else
        weights = dist_weights;
    end

    % Valid indices
    x_valid = x_grid(valid_mask);
    y_valid = y_grid(valid_mask);

    % Update maps
    lin_idx = sub2ind(size(heatmap), y_valid, x_valid);
    heatmap(lin_idx) = heatmap(lin_idx) + weights * score;
    count_map(lin_idx) = count_map(lin_idx) + weights;
end
