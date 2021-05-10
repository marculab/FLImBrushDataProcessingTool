function [out,scale] = processImg(im, posData, FLImData, radius, alpha, SNR_low)

% settings
rad = radius; %15; %3; %10; % radius
% SNR filters
snr_lowerbound = SNR_low;
% color scale
scale = cell(3,1);
out = cell(3,1);
for i =1:3
low = floor(min(FLImData.lt{i})/0.5)*0.5;
high = ceil(max(FLImData.lt{i}/0.5)*0.5);
scale{i} = [low high];
end

for dest_channel = 1:3
fprintf('Plot Lifetime Channel %d\n', dest_channel);
fprintf('SNR LowerBound = %f, Radius = %f, Alpha = %f', snr_lowerbound, rad, alpha);

% setup video
%set(gcf,'Visible', 'off');
jet_cmap =  jet;
%ssx = vr.Width;
%ssy = vr.Height;
ssx = size(im, 2);
ssy = size(im, 1);
% preallocation
overlay = cell(3,1);
val_field = cell(3,1);
accum = cell(3,1);

for i=1:3
    overlay{i}  = uint8(zeros(ssy,ssx,3));
    val_field{i} = double(zeros(ssy,ssx,1));
    accum{i} = double(zeros(ssy,ssx,1));
end

%df1 = im;
message = sprintf('Processing channel %d',dest_channel);
f = waitbar(0,message); % creat waitbar
numOfFrames = length(posData.frames);
BarStep = round(0.01*numOfFrames);
for i = 1:numOfFrames
    % get current video frame
    %df1= vr.read(i);
    % get current segmentation position
    px = posData.px(i);
    py = posData.py(i);
    
    current_value = FLImData.lt{dest_channel}(i);
    current_snr = FLImData.snr{dest_channel}(i);
    % conditions to skip
    if px == 0 || py == 0 || isnan(current_value) || current_value == 0
        continue
    end
    if isnan(current_snr) || current_snr < snr_lowerbound
        continue;
    end
    scaleLow = scale{dest_channel}(1);
    scaleHigh = scale{dest_channel}(2);
    if current_value<scaleLow
        current_value = scaleLow;
    end
    if current_value>scaleHigh
        current_value = scaleHigh;
    end
    ind1 = ceil((current_value-scaleLow)/(scaleHigh-scaleLow)*256);
    [overlay{dest_channel}, val_field{dest_channel}, accum{dest_channel}] = drawCirc( [px,py], rad*0.7, rad, overlay{dest_channel}, jet_cmap(ind1,:)*254+1, val_field{dest_channel}, ind1, accum{dest_channel} );
    
    if mod(i,BarStep)
        waitbar(i/numOfFrames,f,sprintf('Processing channel %d %d %%',dest_channel, round(i/numOfFrames*100)));
    end
end
delete(f);
df1 = im;
df1(~(overlay{dest_channel}(:,:,:) == 0)) = alpha*overlay{dest_channel}(~(overlay{dest_channel}(:,:,:) == 0) ) + (1-alpha)*df1( ~(overlay{dest_channel}(:,:,:) == 0) );
out{dest_channel} = df1;
%set(gcf,'Visible', 'off');



% close window
%close(h);
end

