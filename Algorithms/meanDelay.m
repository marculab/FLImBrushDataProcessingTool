function MD = meanDelay(waveform, dt)
% waveform is colume vector
t = 0:size(waveform,1)-1;
MD = (t+0.5)*waveform./sum(waveform);
MD = MD.*dt;
end
