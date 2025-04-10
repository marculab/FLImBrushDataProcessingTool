function [A, T, avglife, intensity, fitt, raw, decay] = multiexp_fit(spec,dt,laser,orders,taus)
% spec is waveform
% dt is time resolution
% laser is irf
% orders is number of exponentials
% taus are fixed tau values, leave empty '[]' if using auto tau, range of tau is [0.01 25] ns

if ~isempty(taus) % if taus are specified
    if length(taus)~=orders % if taus do not match orders, throw error
        error('Size of taus does not match orders')
    end
end
    
switch orders
    case 1
        [A1, T1, avglife, intensity, fitt, raw, decay] = monoexp_fit(spec,dt,laser,taus);
        A = A1;
        T = T1;
    case 2
        [A1, A2, T1, T2, avglife, intensity, fitt, raw, decay] = biexp_fit(spec,dt,laser,taus);
        A = [A1' A2']';
        T = [T1' T2']';
    case 3
        [A1, A2, A3, T1, T2, T3, avglife, intensity, fitt, raw, decay] = triexp_fit(spec,dt,laser,taus);
        A = [A1' A2' A3']';
        T = [T1' T2' T3']';
    case 4
        [A1, A2, A3, A4, T1, T2, T3, T4, avglife, intensity, fitt, raw, decay] = quadriexp_fit(spec,dt,laser,taus);
        A = [A1' A2' A3' A4']';
        T = [T1' T2' T3' T4']';
    otherwise
        ws = sprintf('Multi-exponential fit with %d orders was not implemented!',orders);
        warning(ws)
end

end