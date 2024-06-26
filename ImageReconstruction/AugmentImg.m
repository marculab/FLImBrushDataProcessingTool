function [augmentedImg,scale] = AugmentImg(im, posData, dataToAugment, scale, radius, alpha, colormap_in)

% fprintf('Radius = %f, Alpha = %f\n', radius, alpha);

% setup video
%set(gcf,'Visible', 'off');
C_map =  colormap_in;
%ssx = vr.Width;
%ssy = vr.Height;
ssx = size(im, 2);
ssy = size(im, 1);
% preallocation
overlay  = uint8(zeros(ssy,ssx,3));
val_field = double(zeros(ssy,ssx,1));
accum = double(zeros(ssy,ssx,1));

%df1 = im;

f = waitbar(0,'Starting'); % creat waitbar
numOfFrames = length(posData.px);
BarStep = round(0.01*numOfFrames);
for i = 1:numOfFrames
    % get current video frame
    %df1= vr.read(i);
    % get current segmentation position
    px = posData.px(i);
    py = posData.py(i);
    
    current_value = dataToAugment(i);
    
    % conditions to skip
    if px == 0 || py == 0 || isnan(current_value) || current_value == 0
        continue
    end
    scaleLow = scale(1);
    scaleHigh = scale(2);
    if current_value<scaleLow
        current_value = scaleLow;
    end
    if current_value>scaleHigh
        current_value = scaleHigh;
    end
    ind1 = ceil((current_value-scaleLow)/(scaleHigh-scaleLow)*255)+1; % make sure index is between 1 to 256
    [overlay, val_field, accum] = drawCirc( [px,py], radius*0.7, overlay, ind1, val_field, accum, C_map);
    
    if mod(i,BarStep)
        waitbar(i/numOfFrames,f,sprintf('Processing %d %%', round(i/numOfFrames*100)));
    end
end
    delete(f);
    augmentedImg = im;  
    mask2D = ~(sum(overlay,3) == 0);
    mask = cat(3,mask2D,mask2D,mask2D); 
    augmentedImg(mask) = alpha*overlay(mask) + (1-alpha)*augmentedImg(mask); 
end

