function [ F ] = processVideo( vr, posData, ltData, filename)

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
if ~plot_redox
    fprintf('Plot Lifetime Channel %d\n', dest_channel);
else
    fprintf('Plot Redox Ratio\n');
end
fprintf('SNR LowerBound = %f, Radius = %f, Alpha = %f', snr_lowerbound, rad, alpha);

% setup video
jet_cmap =  jet(256);
ssx = vr.Width;
ssy = vr.Height;
for i=1:3
    overlay{i}  = uint8( zeros(ssy,ssx,3));
    val_field{i} = double( zeros(ssy,ssx,1));
    accum{i} = double( zeros(ssy,ssx,1));
end

% setup figure
h=figure('Position', [100 100 640 300],'Color',[1 1 1]);


c=1;

if ~plot_redox
    vw = VideoWriter([filename 'out_ch', num2str(dest_channel),'.mp4'],'MPEG-4');
else
    vw = VideoWriter(['out_rr','.mp4'],'MPEG-4');
end

open(vw);
OldFrameNum = 0;
for j = 1: max(posData.frames)
    % get current video frame
    currentFrameNum = j;
    if ~(OldFrameNum==currentFrameNum)
        df1= vr.read(currentFrameNum);
        OldFrameNum = currentFrameNum;
    end
    idx = (posData.frames==currentFrameNum);
    idx = find(idx);
    
    if isempty(idx)
        continue
    end
    for k = 1:length(idx)
        i = idx(k);
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
            ind1 = round((current_value-scale_from(dest_channel))/(scale_to(dest_channel)-scale_from(dest_channel))*255+1);
            [overlay{dest_channel}, val_field{dest_channel}, accum{dest_channel}] = drawCirc( [px,py], rad*0.7, rad, overlay{dest_channel}, jet_cmap(ind1,:)*254+1, val_field{dest_channel}, ind1, accum{dest_channel} );
            df1( ~(overlay{dest_channel}(:,:,:) == 0) ) = alpha*overlay{dest_channel}( ~(overlay{dest_channel}(:,:,:) == 0) ) + (1-alpha)*df1( ~(overlay{dest_channel}(:,:,:) == 0) );
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
            ind1 = round((current_redox-scale_redox_from)/(scale_redox_to-scale_redox_from)*63+1);
            % borrow overlay{1} to store overlaid redox ratio
            [overlay{1}, val_field{1}, accum{1}] = drawCirc( [px,py], rad*0.7, rad, overlay{1}, jet_cmap(ind1,:)*254+1, val_field{1}, ind1, accum{1} );
            df1( ~(overlay{1}(:,:,:) == 0) ) = alpha*overlay{1}( ~(overlay{1}(:,:,:) == 0) ) + (1-alpha)*df1( ~(overlay{1}(:,:,:) == 0) );
        end
    end
    
    %set(gcf,'Visible', 'off');
    if mod(j,max(posData.frames))==0
        
        imshow(df1);
        title(sprintf('Frame %d',currentFrameNum))
        colormap(h, jet);
        
        if ~plot_redox
            caxis([scale_from(dest_channel) scale_to(dest_channel)]);
            h0 = colorbar;
            ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
        else
            caxis([scale_redox_from scale_redox_to]);
            h0 = colorbar;
            ylabel(h0, 'Redox Ratio')
        end
        
        set(gca,'LooseInset',get(gca,'TightInset'))
        
        % save overlaid frame for video export
        %disp(c);
        f = getframe(gcf);
        fprintf('Frame %d\n ', j)
        writeVideo(vw, f);
    end
    c=c+1;
end

% close window
close(h);
close(vw)

end

