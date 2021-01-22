function [G,S] = PhasorTransform(data, w, dt)

G = 0;
S = 0;
area = 0;
for bin = 1:(length(data) - 1)
    
    G = G + data(bin).*cos(w*dt.*(bin-.5))*dt;
    S = S + data(bin).*sin(w*dt.*(bin-.5))*dt;
    area = area + (data(bin) + data(bin+1)).*dt./2;
    
end

G = G./area;
S = S./area;

end