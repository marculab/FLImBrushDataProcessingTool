% code to interplocation of V6 data
%%
close all
clear all
clc

addpath(genpath(pwd))

%% load in data
root = 'C:\Users\Xiangnan\Box Sync\V6 Clinical Data\20210505_neuro024';
[DeConfile,DeConpath] = uigetfile([root '\*.mat'],'Please select DeCon file');
[Txtfile,Txtpath] = uigetfile([root '\*.txt'],'Please select text file');
[Posfile,Pospath] = uigetfile([root '\*.mat'],'Please select CNN Segementation Pos file');
[videoName,videopath] = uigetfile([root '\*.avi'],'Please select video file');

%% load in data
load(fullfile(DeConpath,DeConfile))
load(fullfile(Pospath,Posfile))
VideoData = importVideoTextFile(fullfile(Txtpath,Txtfile));
RepRate = Ch1DataObj.laserRepRate;
%% calculate num of data points
% shift = length(Ch1LT)-(MetaData(end,9)-MetaData(1,9))/1000*120/4;
% shift = round(shift);
shift=0; %shift need to be 0!

Ch1INTCorr = circshift(Ch1DataObj.INTsAllGainCorrected,shift);
Ch1LT = circshift(Ch1DataObj.LTsAll,shift);
Ch1SNR = circshift(Ch1DataObj.SNR,shift);
G1 = circshift(Ch1DataObj.gain,shift);

Ch2INTCorr = circshift(Ch2DataObj.INTsAllGainCorrected,shift);
Ch2LT = circshift(Ch2DataObj.LTsAll,shift);
Ch2SNR = circshift(Ch2DataObj.SNR,shift);
G2 = circshift(Ch2DataObj.gain,shift);

Ch3INTCorr = circshift(Ch3DataObj.INTsAllGainCorrected,shift);
Ch3LT = circshift(Ch3DataObj.LTsAll,shift);
Ch3SNR = circshift(Ch3DataObj.SNR,shift);
G3 = circshift(Ch3DataObj.gain,shift);

%% repopulate output data mat with interplation
% reset start time to 0
VideoData(:,9) = VideoData(:,9)-VideoData(1,9);
% use CNN segmentation location data
try
    VideoData(:,6:7) = double(pos(:,7:8));
catch
    VideoData(:,6:7) = double(pos_new(:,1:2));
end
% Length = 70000;
Length = VideoData(end,9);
VideoData(VideoData(:,9)>Length,:) = [];
numOfVideoFrame = size(VideoData,1);
removeFlag = zeros(numOfVideoFrame,1);
%------------remove single 0 lines---------------------------------------%
for i = 1:numOfVideoFrame
    posCurrent = VideoData(i,6)+ VideoData(i,7);
    if posCurrent==0 % if current pos is all 0 check before and after
        if i==1
            posBefore = 1;
            posAfter = VideoData(i+1,6)+ VideoData(i+1,7);
        else if i==numOfVideoFrame
                posBefore = VideoData(i-1,6)+ VideoData(i-1,7);
                posAfter = 1;
            else
                posBefore = VideoData(i-1,6)+ VideoData(i-1,7);
                posAfter = VideoData(i+1,6)+ VideoData(i+1,7);
            end
            removeFlag(i) = posBefore&posAfter;
        end
    end
end
removeIdx = find(removeFlag);
VideoData(removeIdx,:) = [];

% remove all 0 points
temp = VideoData(:,6)+VideoData(:,7);
ZeroIdx = find(temp==0);
VideoData(ZeroIdx,:)=[];

frameNum = 1:size(VideoData,1);
frameNum = frameNum';
outputMat = zeros(size(Ch1LT,1),23);

time = 0:1:size(Ch1LT,1)-1;
time = time*1000/RepRate*4;
frameT = VideoData(:,9);
[frameT,ia] = unique(frameT); % find duplicate frame
VideoData = VideoData(ia,:); % remove duplicated data
frameNum = frameNum(ia); % remove duplicated frame
frameIdx = interp1(frameT,frameNum,time);
frameIdx = ceil(frameIdx)';
outputMat(:,1) = frameIdx;
VideoData = [frameNum VideoData];


%% interpolate locations
xx = interp1(VideoData(:,10),VideoData(:,7),time);
xx(isnan(xx))=0; % replace NaN with 0
yy = interp1(VideoData(:,10),VideoData(:,8),time);
yy(isnan(yy))=0; % replace NaN with 0
rr = interp1(VideoData(:,10),VideoData(:,9),time);
rr(isnan(rr)) = 0;
outputMat(:,7) = xx;
outputMat(:,8) = yy;
outputMat(:,9) = rr;
outputMat(:,12) = Ch1LT;
outputMat(:,16) = Ch1INTCorr;
outputMat(:,13) = Ch2LT;
outputMat(:,17) = Ch2INTCorr;
outputMat(:,14) = Ch3LT;
outputMat(:,18) = Ch3INTCorr;
outputMat(:,20) = Ch1SNR;
outputMat(:,21) = Ch2SNR;
outputMat(:,22) = Ch3SNR;
outputMat((xx+yy)==0,:)=[]; % remove bad data points

%% save data
cd(DeConpath)
[filepath,name,ext] = fileparts(DeConfile);
save([name '_ImgRecon.mat'],'outputMat')
disp('Finished')
close all

%%
cd(DeConpath)
% replotVideo(['videos\' videoName '.avi'], [name '_interp.mat'])
replotInterpImage([name '_ImgRecon.mat'],fullfile(videopath,videoName),10,0.3,20)
