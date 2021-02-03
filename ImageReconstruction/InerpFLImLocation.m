% code to interplocation of V6 data
%%
close all
clear all
clc

addpath(genpath(pwd))

%% load in data
root = 'C:\Users\Xiangnan\Documents\MyGitRepo\FLImBrushDataProcessingTool\TestData\20200921Lamp1mW120Hz';
[DeConfile,DeConpath] = uigetfile([root '\*.mat'],'Please select DeCon file');
[Posfile,Pospath] = uigetfile([root '\*.mat'],'Please select CNN Segementation Pos file');
[Txtfile,Txtpath] = uigetfile([root '\*.txt'],'Please select text file');
[videoName,videopath] = uigetfile([root '\*.avi'],'Please select video file');
RepRate = 120;

%% load in data
load(fullfile(DeConpath,DeConfile))
load(fullfile(Pospath,Posfile))
MetaData = importVideoTextFile(fullfile(Txtpath,Txtfile));

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
MetaData(:,9) = MetaData(:,9)-MetaData(1,9);
% Length = 70000;
Length = MetaData(end,9);
MetaData(MetaData(:,9)>Length,:) = [];
numOfData = size(MetaData,1);
removeFlag = zeros(numOfData,1);
%------------remove single 0 lines---------------------------------------%
for i = 1:numOfData
    posCurrent = MetaData(i,6)+ MetaData(i,7);
    if posCurrent==0 % if current pos is all 0 check before and after
        if i==1
            posBefore = 1;
            posAfter = MetaData(i+1,6)+ MetaData(i+1,7);
        else if i==numOfData
                posBefore = MetaData(i-1,6)+ MetaData(i-1,7);
                posAfter = 1;
            else
                posBefore = MetaData(i-1,6)+ MetaData(i-1,7);
                posAfter = MetaData(i+1,6)+ MetaData(i+1,7);
            end
            removeFlag(i) = posBefore&posAfter;
        end
    end
end
        removeIdx = find(removeFlag);
        MetaData(removeIdx,:) = [];
        
        % replace all 0 points with NaN
        Idx = MetaData(:,6)+MetaData(:,7);
        temp = MetaData(:,6:7);
        temp(temp==0)=NaN;
        MetaData(:,6:7)=temp;
        
        frameNum = 1:size(MetaData,1);
        frameNum = frameNum';
        outputMat = zeros(size(Ch1LT,1),23);
        
        time = 0:1:size(Ch1LT,1)-1;
        time = time*1000/RepRate*4;
        frameT = MetaData(:,9);
        [frameT,ia] = unique(frameT); % find duplicate frame
        MetaData = MetaData(ia,:); % remove duplicated data
        frameNum = frameNum(ia); % remove duplicated frame
        frameIdx = interp1(frameT,frameNum,time);
        frameIdx = ceil(frameIdx)';
        outputMat(:,1) = frameIdx;
        MetaData = [frameNum MetaData];
        
        
        %%
        xx = interp1(MetaData(:,10),MetaData(:,7),time);
        xx(isnan(xx))=0; % replace NaN with 0
        yy = interp1(MetaData(:,10),MetaData(:,8),time);
        yy(isnan(yy))=0; % replace NaN with 0
        rr = interp1(MetaData(:,10),MetaData(:,9),time);
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
        cd(DeConpath)
        %% save data
        [filepath,name,ext] = fileparts(DeConfile);
        save([name '_ImgRecon.mat'],'outputMat')
        disp('Finished')
        close all
        
        %%
        cd(DeConpath)
        
        % replotVideo(['videos\' videoName '.avi'], [name '_interp.mat'])
        replotInterpImage(videoName, [name '_ImgRecon.mat'])
