function [PhasorS]=ProcessPhasor(Data,Channels,Harmonics)

Sizecheck= size(Data{1,1});
if (Sizecheck(1)== 0 )
PhasorS=zeros(7,0);
    return
end
% If Data is empty, it returns 

%Gets Data structire and returns phasors in cell 
[DimT,DimX]=size(Data{1,1}.channeldata.data);



PhasorS=zeros(7,DimX);


ir=ProcessIr(Data);

if (Channels==4) 
  irE=cat(1,ir(:,1),ir(:,2),ir(:,3),ir(:,4));
else
  irE=cat(1,ir(:,1),ir(:,2),ir(:,3));
end

decay=zeros(DimT,Channels);


for x=1:DimX
    for Channel=1:Channels
    Z=(Data{Channel,1}.channeldata.data(:,x));
    %Z=ShiftArray(Z,Data{Channel,1}.shift);
    %Cropping zeros
    %Z=Z.*(Z>0);
%     
%      BG=sum(Z(580:680,:))/100;
%      Z=Z+BG;
    decay(:,Channel)=Z;
    PhasorS(Channel,x)=ComputePhasor(decay(:,Channel),ir(:,Channel),Harmonics,1);
    end
    
%Creating extended decay
if (Channels==4)
decayE=cat(1,decay(:,1),decay(:,2),decay(:,3),decay(:,4));
Spectrum=zeros(4,1);
Spectrum(1)=sum(decay(:,1));
Spectrum(2)=sum(decay(:,2));
Spectrum(3)=sum(decay(:,3));
Spectrum(4)=sum(decay(:,4));
Spectrum=Spectrum/sum(Spectrum);
else
decayE=cat(1,decay(:,1),decay(:,2),decay(:,3));
Spectrum=zeros(3,1);
Spectrum(1)=sum(decay(:,1));
Spectrum(2)=sum(decay(:,2));
Spectrum(3)=sum(decay(:,3));
Spectrum=Spectrum/sum(Spectrum);

end
%ST Phasor
irE(:)=0;
PhasorS(5,x)=ComputePhasor(decayE,irE,Harmonics,0);

DecayF=zeros(Channels,1);
for i = 1 : Channels
DecayF(i)= PhasorS(Channel,x)*Spectrum(i);
end
PhasorS(6,x)=(ComputePhasor(DecayF,irE,1,0)+(1+1i))/2;

% Spectral Phasor
PhasorS(7,x)=(ComputePhasor(Spectrum,irE,1,0)+(1+1i))/2;

end
end
