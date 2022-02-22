function [y]=ComputePhasor(S,Ir,Harmonics,Normalize)

Harmonics=Harmonics+1;

Decay=fft(S);

if (sum(Ir)~=0)
    %Ir=Ir.*(Ir>0);
IR=fft(Ir);
IR(Harmonics)=IR(Harmonics)./IR(1);
else
IR(Harmonics)=1;    
end


if (Normalize==1)
Decay(Harmonics)=Decay(Harmonics)./(Decay(1));
end 
y=conj(Decay(Harmonics)/IR(Harmonics));
end
