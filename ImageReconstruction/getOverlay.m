function [overlay,overlayAllFrame] = getOverlay(imgSize, posData, dataToAugment, scale, radius)
% generate overlay data only
% fprintf('Radius = %f', radius);
% overlay is the final render
% overlayAllFrame contains the rendering for each video frame

frameIdx = posData.frameIdx;
frameIdx(isnan(frameIdx)) = []; % remove NaN
frameIdx = [frameIdx;frameIdx(end)+1]; % append one data to makes sure last frame is included

frameIdxDiff = diff(frameIdx);
frameChangeIdx = find(frameIdxDiff==1);
overlayAllFrame = cell(size(frameChangeIdx)); % cell array to store all overlay
% setup video
%set(gcf,'Visible', 'off');
jet_cmap =  jet(256); % colormap
%ssx = vr.Width;
%ssy = vr.Height;
ssx = imgSize(2);
ssy = imgSize(1);
% preallocation
overlay  = uint8(zeros(ssy,ssx,3)); % initilize
val_field = double(zeros(ssy,ssx,1)); % initilize
accum = double(zeros(ssy,ssx,1)); % initilize

%df1 = im;
numOfPoints = length(posData.px);
scaleLow = scale(1);
scaleHigh = scale(2);
% f = waitbar(0,'Starting'); % creat waitbar
BarStep = round(0.01*numOfPoints);
for i = 1:numOfPoints
    i
    % get current video frame
    %df1= vr.read(i);
    % get current segmentation position
    px = posData.px(i);
    py = posData.py(i);
    
    current_value = dataToAugment(i);
    
    % conditions to skip
    if isnan(px) || isnan(py) || isnan(current_value) || current_value == 0
        % save old overlay to cell if current data is not valid
        if sum(i==frameChangeIdx)
            frameTemp = frameIdx(i);
            overlayAllFrame{frameTemp} = overlay;
        end
        continue
    end
    % if not skip
    if current_value<scaleLow
        current_value = scaleLow;
    end
    if current_value>scaleHigh
        current_value = scaleHigh;
    end
    ind1 = ceil((current_value-scaleLow)/(scaleHigh-scaleLow)*255)+1; % find colormap index, range 1 to 256
    [overlay, val_field, accum] = drawCirc( [px,py], radius*0.7, overlay, ind1, val_field, accum, jet_cmap); % update overlay
    % if last points for the frame save overlay to cell
    if sum(i==frameChangeIdx)
        frameTemp = frameIdx(i);
        overlayAllFrame{frameTemp} = overlay;
    end
    % update progress bar
    % if mod(i,BarStep)
    %     waitbar(i/numOfPoints,f,sprintf('Processing %d %%', round(i/numOfPoints*100)));
    % end
end
% delete(f);
end

