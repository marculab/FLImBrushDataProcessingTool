function saveDeconLite(Ch1DataObjIn,Ch2DataObjIn,Ch3DataObjIn,Ch4DataObjIn,saveFileFullPath)
if isempty(Ch4DataObjIn)
    numOfChannel = 3;
else
    numOfChannel = 4;
end
Ch1DataObj = struct;
Ch2DataObj = struct;
Ch3DataObj = struct;
Ch4DataObj = struct;

VarListLite={'laserRepRate','alpha', 'dataT', 'gain','irfIdx','Lg_INTs',...
    'Lg_INTsGainCorrected', 'K', 'LaguerreBasis', 'Lg_LCs','Lg_LTs',...
    'noise','SNR','Chi2','Ph_H1S','Ph_H1G','Ph_H2S','Ph_H2G',...
    'Ph_H3S','Ph_H3G','Ph_H4S','Ph_H4G'};

VarListDecon={'laserRepRate','alpha', 'dataT', 'gain','irfIdx','Lg_INTs',...
    'Lg_INTsGainCorrected', 'K', 'LaguerreBasis', 'Lg_LCs','Lg_LTs',...
    'noise','SNR','stat_test.chi2.stat''','Ph_H1S','Ph_H1G','Ph_H2S','Ph_H2G',...
    'Ph_H3S','Ph_H3G','Ph_H4S','Ph_H4G'};

for ii=1:numOfChannel
    for jj=1:length(VarListLite)
        eval(['Ch' num2str(ii) 'DataObj.' VarListLite{jj} '=Ch' num2str(ii) 'DataObjIn.' VarListDecon{jj} ';'])
        % eval(['SP_G' '=data.SP_G'';'])
        % eval(['SP_S' '=data.SP_S'';'])
    end
end
[filepath,name,ext] = fileparts(saveFileFullPath);
name = strcat(name,'_lite',ext);
new_filePath = fullfile(filepath,name);

if numOfChannel == 3
    save(new_filePath,'Ch1DataObj','Ch2DataObj','Ch3DataObj')
else
    save(new_filePath,'Ch1DataObj','Ch2DataObj','Ch3DataObj','Ch4DataObj')
end