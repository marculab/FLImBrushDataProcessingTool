function [Gout,Sout] = MPcorrection(G,S,M,P)
Gout = (G.*cos(P) + S.*sin(P))./M;
Sout = (-G.*sin(P) + S.*cos(P))./M;
end