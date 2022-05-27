classdef ChannelData
    % CHANNELDATA class, store raw data and relavent parameters for each spectral channel.
    properties
        data % raw data from one spectral channel, DC and BG removed and truncated, 2D matrix
        iIRF % instrument response funciton of spectral channel, 1D vector
        dt % time resolution of ddigitizer, scaler
        bw % bandwidth of the system, scaler
        decon_idx % index of waveforms that is being processed, after saturation and low signal removal
        noise % noise of the waveform, calculated using the last 100 points of raw wavefrom, before trunction
        SNR % singal to noise ratio of each wavefrom, 1D vector
        gain_list % 1D vector of gain of each wavefrom
    end
    methods
        function obj = ChannelData(data,iIRF,dt,bw,deconidx,noise,gain_list)
            % constructor, create object to store waveforms and parameters from one spectral channel
            obj.data = data;
            obj.iIRF = iIRF;
            obj.dt = dt;
            obj.bw = bw;
            obj.decon_idx = deconidx;           
            obj.gain_list = gain_list;
            % calculate SNR
            if ~isempty(noise)
                obj.noise = noise;
                data_max = max(data,[],1);
                obj.SNR = 20*log10(data_max./noise);
                obj.SNR(data_max < noise) = 0;
                obj.SNR = obj.SNR';
            else
                obj.noise = NaN;
            end
        end
    end
    
    
    
end