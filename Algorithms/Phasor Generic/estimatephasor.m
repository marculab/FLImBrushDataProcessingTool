function [G, S] = estimatephasor(decay, ref, ref_tau, order, delta_t)
% decay is the array of the decay curve
% ref is the decay curve of the reference or IRF curve
% ref_tau is the lifetime of the reference compound (zero for IRF)
% freq is the freq domain freq you want to use to plot the phasor
% delta_t is the size (time) of each of the bins of the decay curve

w = order * 2 * pi / (size(decay,1) * delta_t);

% calculate ref phasor
[G_ref, S_ref] = PhasorTransform(ref,w,delta_t);


% calculate phase and modulation corrections
[M_cor,ph_cor] = GS2MP(G_ref,S_ref,w*ref_tau);


% calculate data phasor
G = zeros(size(decay,2),1);
S = zeros(size(decay,2),1);
for i = 1: size(decay,2)
    [Gdec, Sdec] = PhasorTransform(decay(:,i),w,delta_t);
    [G(i),S(i)] = MPcorrection(Gdec,Sdec,M_cor,ph_cor);
end
end