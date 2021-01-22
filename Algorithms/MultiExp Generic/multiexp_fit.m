function [A, T, avglife, intensity, fitt, raw] = multiexp_fit(spec,dt,laser,orders)
switch orders
    case 1
        [A1, T1, avglife, intensity, fitt, raw] = monoexp_fit(spec,dt,laser);
        A = A1;
        T = T1;
    case 2
        [A1, A2, T1, T2, avglife, intensity, fitt, raw] = biexp_fit(spec,dt,laser);
        A = [A1' A2']';
        T = [T1' T2']';
    case 3
        [A1, A2, A3, T1, T2, T3, avglife, intensity, fitt, raw] = triexp_fit(spec,dt,laser);
        A = [A1' A2' A3']';
        T = [T1' T2' T3']';
    case 4
        [A1, A2, A3, A4, T1, T2, T3, T4, avglife, intensity, fitt, raw] = quadriexp_fit(spec,dt,laser);
        A = [A1' A2' A3' A4']';
        T = [T1' T2' T3' T4']';
    otherwise
        ws = sprintf('Multi-exponential fit with %d orders was not implemented!',orders);
        warning(ws)
end

end