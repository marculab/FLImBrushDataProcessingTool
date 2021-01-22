function [out_BG, BG_filename]= loadTRIPLEXBG(BGRootPath)

[filename_BG, foldername_BG] = uigetfile('*.*', 'Please Select the Background Data from TRIPLEX', BGRootPath);
if ~foldername_BG
    
    out_BG = [];
    return
    
end
% try loading BG in different format
try
    BG = load(fullfile(foldername_BG, filename_BG))';
    
catch err0
    
    try
        
        BG = loadBinary(fullfile(foldername_BG, filename_BG), 'int16', 0)';
        
    catch err1
        
        try
            
            BG = -LoadBackground(fullfile(foldername_BG, filename_BG))';
            
        catch err2
            
            BG = loadBinary(fullfile(foldername_BG, filename_BG), 'double', 0)';
            
        end
        
    end
    
end

ssr = mean(BG((end - 100):end, :));
BG = bsxfun(@minus, BG, ssr);
% BG = bsxfun(@rdivide, BG, abs(min(BG))); %normalise BG as BG is negative
BG = mean(BG, 2);
BG_filename = fullfile(foldername_BG,filename_BG);
out_BG = BG;