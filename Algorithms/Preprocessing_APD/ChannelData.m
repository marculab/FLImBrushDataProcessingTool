classdef ChannelData
    properties
        data
        iIRF
        dt
        bw
        decon_idx
        noise
        SNR
        gain_list
    end
    methods
        function obj = ChannelData(data,iIRF,dt,bw,deconidx,noise,gain_list)
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