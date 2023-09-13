% code to add phasor data to the reporcessed HN Data file
%%
clear
close all
clc
addpath(genpath('..'))
%%  load aggregate data
HNAggregateInputDataRunLevel = importHNRunData('H&N Aggregate Input Data - Run_Level.csv');
HNAggregateInputDataRunLevel = HNAggregateInputDataRunLevel(HNAggregateInputDataRunLevel.Patient>100,:);
numOfRuns = size(HNAggregateInputDataRunLevel,1);
dataRoot = 'Y:\V2 Sacramento Database\Da Vinci Robot Study (100 patients)\Data_100_Patient_Study';
for i = 1:numOfRuns
    runfile = HNAggregateInputDataRunLevel.DeconvolutionFile(i)
    if ~(runfile=="")
        runfilePath = fullfile(dataRoot,runfile);
        variableInfo = who('-file', runfilePath);
        load(runfilePath)
        for H = 1:4 % run phasor for harmonic of 1 to 4
            runPhasor(Ch1DataObj, H);
            runPhasor(Ch2DataObj, H);
            runPhasor(Ch3DataObj, H);
        end

        CH1WFGainCorr = get(Ch1DataObj,'GainCorrectedWF');
        CH2WFGainCorr = get(Ch2DataObj,'GainCorrectedWF');
        CH3WFGainCorr = get(Ch3DataObj,'GainCorrectedWF');
        %-------------------------------------------compute extended phasor-----------------------------------------------------
        WFExtended = [CH1WFGainCorr;CH2WFGainCorr;CH3WFGainCorr];
        numOfWF = size(WFExtended,2);
        EOP = zeros(numOfWF,1);
        for j = 1:numOfWF
            EOP(i) = ComputePhasor(WFExtended(:,i),0,1,0);
        end
        EOP_H1G = real(EOP);
        EOP_H1S = imag(EOP);
        %-------------------------------------------compute spectral phasor--------------------------------------------
        Spectrum = zeros(numOfWF,3);
        Spectrum(:,1) = sum(CH1WFGainCorr)';
        Spectrum(:,2) = sum(CH2WFGainCorr)';
        Spectrum(:,3) = sum(CH3WFGainCorr)';
        SP = zeros(numOfWF,1);
        for j = 1: numOfWF
            SP(i) = (ComputePhasor(Spectrum(i,:),0,1,0)+(1+1i))/2;
        end
        SP_G = real(SP);
        SP_S = imag(SP);

        save(runfilePath, 'dataInfoObj','Ch1DataObj','Ch2DataObj','Ch3DataObj','EOP_H1G','EOP_H1S','SP_G','SP_S','-v7.3');
    end
end