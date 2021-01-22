function [M,P] = GS2MP(G,S,wtau)
M = sqrt(G.^2 + S.^2)./sqrt(1./(1+wtau.^2));
P = atan2(S,G) + atan(wtau);
end