function [ img, val_field, accum ] = drawCirc( c, r, N, img, val, val_field, vls, accum )
%DRAWCIRC Summary of this function goes here
%   Detailed explanation goes here
N=100;
[imgH, imgW,~] = size(img);
% create a circular mask
t = linspace(0, 2*pi, N);
BW = poly2mask(r*cos(t)+c(1), r*sin(t)+c(2), imgH, imgW);

if vls>0
    %val_field(BW) = ((val_field(BW)./max(uint8(BW(BW)),accum(BW)) ) + vls)./(accum(BW)+1);
    val_field(BW) = ((val_field(BW).*accum(BW) ) + vls)./(accum(BW)+1);
    accum(BW) = accum(BW)+1;
end

jet_cmap = jet;
% overlay filled circular shape by using the mask
% to fill the image with the desired color (for all three channels R,G,B)
clr = val; %[0 0 255];            % blue color
a = 1;                    % blending factor
z = false(size(BW));
mask = cat(3,BW,z,z); 
img(mask) = jet_cmap( round(val_field(BW)),1 )*254+1; %a*clr(1)  + (1-a)*img(mask);
mask = cat(3,z,BW,z); 
img(mask) = jet_cmap( round(val_field(BW)),2 )*254+1; %img(mask) = a*clr(2)  + (1-a)*img(mask);
mask = cat(3,z,z,BW); 
img(mask) = jet_cmap( round(val_field(BW)),3 )*254+1; %img(mask) = a*clr(3)  + (1-a)*img(mask);

end

