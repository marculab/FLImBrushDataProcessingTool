function [out_iRF,iRF_name] = loadTRIPLEXiIRF(iRFRootPath)

[filename, foldername] = uigetfile('*.mat;*.dat', 'Please Select the TRIPLEX iIRF File', iRFRootPath);

if ~foldername
    
    out_iRF = [];
    iRF_name = [];
    return;
    
end

full_path = fullfile(foldername, filename);
iIRF_in = load(full_path);
try
fName = fieldnames(iIRF_in);
iIRF = getfield(iIRF_in, char(fName));

% normalize iRF area under curve to 1
for i = 1:length(iIRF)
    iIRF{i} = iIRF{i}./sum(iIRF{i});
end
catch
    iIRF{1} = iIRF_in;
end
out_iRF = iIRF;
iRF_name=filename;
return

