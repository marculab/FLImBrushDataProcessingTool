classdef IRF < handle

properties (Access = protected)
    
end

properties (Access = public)
    ch1; %1D n by 1
    ch2; %1D n by 1
    ch3; %1D n by 1
    ch4; %1D n by 1
    
    ch1Norm;
    ch2Norm;
    ch3Norm;
    ch4Norm;
    
    pathToIRF ; % full path to IRF file
end

methods
   %constructor
   function obj = IRF(irfPath)
       
       obj.pathToIRF = irfPath;
       
       temp = load(irfPath);
       obj.ch1 = temp.laser{1};
       obj.ch2 = temp.laser{2};
       obj.ch3 = temp.laser{3};
       obj.ch4 = temp.laser{4};
       
       obj.ch1Norm = obj.ch1/sum(obj.ch1);
       obj.ch2Norm = obj.ch2/sum(obj.ch2);
       obj.ch3Norm = obj.ch3/sum(obj.ch3);
       obj.ch4Norm = obj.ch4/sum(obj.ch4);
   end
     
end


end