function [] = processImgLocation( im, posData, ltData ,filename)

% settings
rad = 25; %15; %3; %10; % radius
alpha = 0.5; % alpha value for overlay
dest_channel = 2;
% SNR filters
snr_lowerbound = 10;
% color scale
scale_from=[]; scale_to=[];
scale_from(1) = 3.5;
scale_to(1) = 5;
scale_from(2) = 3;
scale_to(2) = 6;
scale_from(3) = 2.5;
scale_to(3) = 5;

% if configuration file exists
if exist('plotConfig.csv', 'file') == 2
    fprintf('plotConfig.csv file exists\n');
    configMat = csvread('plotConfig.csv',1,0);
    dest_channel = configMat(1,1);
    plot_redox = configMat(1,2);
    snr_lowerbound = configMat(1,3);
    rad = configMat(1,4);
    alpha = configMat(1,5);
    scale_from(1) = configMat(1,6);
    scale_to(1) = configMat(1,7);
    scale_from(2) = configMat(1,8);
    scale_to(2) = configMat(1,9);
    scale_from(3) = configMat(1,10);
    scale_to(3) = configMat(1,11);
    scale_redox_from = configMat(1,12);
    scale_redox_to = configMat(1,13);
else
    fprintf('plotConfig.csv file does not exist, using default values\n');
end
rad = 5;
alpha = 0.8;
if ~plot_redox
    fprintf('Plot Lifetime Channel %d\n', dest_channel);
else
    fprintf('Plot Redox Ratio\n');
end
fprintf('SNR LowerBound = %f, Radius = %f, Alpha = %f', snr_lowerbound, rad, alpha);

% setup video
set(gcf,'Visible', 'off');
jet_cmap =  colormap('gray');
%ssx = vr.Width;
%ssy = vr.Height;
ssx = size(im, 2);
ssy = size(im, 1);
for i=1:3
    overlay{i}  = uint8( zeros(ssy,ssx,3)); 
    val_field{i} = double( zeros(ssy,ssx,1)); 
    accum{i} = double( zeros(ssy,ssx,1)); 
end

% setup figure
h=figure; hold on;

set(h, 'Position', [200 200 640 300]);
set(h,'Color',[1 1 1]);

%df1 = im;

for i = 1: length(posData.frames)
    % get current video frame
    %df1= vr.read(i);
    % get current segmentation position
    px = posData.px(i);
    py = posData.py(i);
    
    current_value = ltData.lt{dest_channel}(i);
    current_snr = ltData.snr{dest_channel}(i);
    current_redox = ltData.lt{2}(i) / (ltData.lt{2}(i) + ltData.lt{3}(i)); % rr = (Int 2/(Int 2 + Int 3))
    % conditions to skip
    if px == 0 || py == 0 || isnan(current_value) || current_value == 0 
        continue
    end
    if isnan(current_snr) || current_snr < snr_lowerbound
        continue;
    end
    if ~plot_redox
        if current_value<scale_from(dest_channel)
            current_value = scale_from(dest_channel);
        end
        if current_value>scale_to(dest_channel)
            current_value = scale_to(dest_channel);
        end
        ind1 = 256;
        [overlay{dest_channel}, val_field{dest_channel}, accum{dest_channel}] = drawCircWhite( [px,py], rad, rad, overlay{dest_channel}, jet_cmap(ind1,:)*254+1, val_field{dest_channel}, ind1, accum{dest_channel} );
    else % plot redox ratio
        if isnan(current_redox)
            continue
        end
        if current_redox<scale_redox_from
            current_redox = scale_redox_from;
        end
        if current_redox>scale_redox_to
            current_redox = scale_redox_to;
        end
        ind1 = round((current_redox-scale_redox_from)/(scale_redox_to-scale_redox_from)*254+1);
        % borrow overlay{1} to store overlaid redox ratio
        [overlay{1}, val_field{1}, accum{1}] = drawCircWhite( [px,py], rad, rad, overlay{1}, jet_cmap(ind1,:)*254+1, val_field{1}, ind1, accum{1} );
    end
        fprintf('%i\n ', i);
end

df1 = im;
df1( ~(overlay{dest_channel}(:,:,:) == 0) ) = alpha*overlay{dest_channel}( ~(overlay{dest_channel}(:,:,:) == 0) ) + (1-alpha)*df1( ~(overlay{dest_channel}(:,:,:) == 0) );

%set(gcf,'Visible', 'off');
hold off; imshow(df1);

colormap(h, jet);

% if ~plot_redox
%     caxis([scale_from(dest_channel) scale_to(dest_channel)]);
%     h0 = colorbar;
%     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
% else
%     caxis([scale_redox_from scale_redox_to]);
%     h0 = colorbar;
%     ylabel(h0, 'Redox Ratio')
% end

set(gca,'LooseInset',get(gca,'TightInset'))

if ~plot_redox
    saveas(gcf, [filename 'outimg_loc','.jpg']);
else
    saveas(gcf, ['outimg_rr','.jpg']);
end

% close window
%close(h);
end

