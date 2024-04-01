function saveDeconLite(Ch1DataObjIn,Ch2DataObjIn,Ch3DataObjIn,saveFileFullPath)
Ch1DataObj = struct;
Ch2DataObj = struct;
Ch3DataObj = struct;
VarList={'laserRepRate','alpha', 'dataT', 'gain','irfIdx','Lg_INTs','Lg_INTsGainCorrected', 'K', 'LaguerreBasis', 'Lg_LCs','Lg_LTs','noise','SNR','Ph_H1S','Ph_H1G','Ph_H2S','Ph_H2G','Ph_H3S','Ph_H3G','Ph_H4S','Ph_H4G'};
for ii=1:3
    for jj=1:length(VarList)
        eval(['Ch' num2str(ii) 'DataObj.' VarList{jj} '=Ch' num2str(ii) 'DataObjIn.' VarList{jj} ';'])
        % eval(['SP_G' '=data.SP_G'';'])
        % eval(['SP_S' '=data.SP_S'';'])
    end
end
[filepath,name,ext] = fileparts(saveFileFullPath);
name = strcat(name,'_lite',ext);
new_filePath = fullfile(filepath,name);
save(new_filePath,'Ch1DataObj','Ch2DataObj','Ch3DataObj')
