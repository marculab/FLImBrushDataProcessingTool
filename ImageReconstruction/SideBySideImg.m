% code to generate side by side images of steak
%%
clear all
close all
clc

%% load in images
WLImg = imread('..\Data\20200921Lamp1mW120Hz\Nearest Neighbor\20200921Lamp1mW120Hz_WL.jpg');
LowRepRateImg = imread('..\Data\20200921Lamp1mW120Hz\Nearest Neighbor\CH2 LT\20200921Lamp1mW120Hz12_interp.png');
HighRepRateImg = imread('..\Data\20200921Lamp0.25mW480Hz\Nearest Neighbor\CH2 LT\20200921Lamp0.25mW480Hz12_interp.png');
LowRepIntermImg = imread('..\Data\20200921Lamp1mW120Hz\videos\20200921Lamp1mW120Hz_run12.avioutimg_ch2.jpg');
HighRepIntermImg = imread('..\Data\20200921Lamp0.25mW480Hz\videos\20200921Lamp0.25mW480Hz_run12.avioutimg_ch2.jpg');


%% 
figure('Position', [200 200 800 450])
tiledlayout(2,3)
nexttile
imshow(WLImg)
title('White Light Image')
nexttile
imshow(LowRepIntermImg)
title('120Hz 1mW beofre nearest neighber interp')
nexttile
imshow(HighRepIntermImg)
title('480Hz 0.24mW beofre nearest neighber interp')

nexttile
% do nothing
nexttile
imshow(LowRepRateImg)
title('120Hz 1mW after nearest neighber interp')
nexttile
imshow(HighRepRateImg)
title('480Hz 0.24mW after nearest neighber interp')