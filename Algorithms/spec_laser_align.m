function out=spec_laser_align(spec,laser,K,n_align_fst,shif_range,shif_c,timeout)
% K - order used for alignment
% n_align_fst - number of pulses used for alignment (with highest intensity)
% shif_range - window for shifting laser
% shif_c - initial guess of laser shift

% Recursive timeout == 60 seconds
if isempty(timeout)
    timeout = clock;
end
% show alignment result flag
show_align_result_flag = 0;
force_non_neg_flag = 0;

M=size(spec,1);
% K=8;
alpha=alpha_up(M,K);

n_align=min(n_align_fst,size(spec,2));%first n_align pulses for alignment

L=Laguerre(M,K,alpha);
% L=[L,ones(length(L),1),];%bkgd fitting
[max_spec,max_spec_I]=sort(max(spec),'descend'); % sort WF by peak value
nstep = round(size(spec,2)/n_align);%calculate steps
spec_t=spec(:,max_spec_I(1:nstep:end));%temporary spec matrix for alignment, evenly pick n_align number of WF

[~, I_max] = max(spec_t(:,1));% find data maximum

if force_non_neg_flag % force non-negative
    temp = spec_t(1:I_max,:);
    temp(temp<0) = 0;
    spec_t(1:I_max,:) = temp;
end

%====find the shift range=============================
if ~exist('shif_range','var')||isempty(shif_range)
    shif_range=20;
end

if ~exist('shif_c','var')||isempty(shif_c) %find initial shift based on peak value
    [~,I_laser]=min(laser<max(laser)*0.5);
    [~,I_spec]=min(spec_t(:,1)<max_spec(1)*0.5);
    shif_c=I_spec-I_laser;
end

shif=shif_c-shif_range/2:shif_c+shif_range/2;
%====================================================

err_shif=zeros(1,length(shif));
% shift iRF, fit and calculate error
for si=1:length(shif)
    %     si %print process
    laser_temp=circshift(laser,[shif(si),0]);
    vv=filter(laser_temp,1,L);
    %     cc=vv\spec_t;
    %     h=L*cc;
    %     spec_t_fitted=vv*cc;
    %     error_temp = spec_t-spec_t_fitted;
    %     err_shif(si)=norm(error_temp,'fro');
    D_mat=conv2(eye(size(spec_t,1)),[1,-3,3,-1],'valid')'*L;
    D=D_mat;
    H=vv'*vv; %positive definite matrix
    H_chol=chol(inv(H)); %Cholesky decomposition
    C=H_chol*D';
    l1=H_chol*vv';
    lam=zeros(size(D,1),size(spec_t,2));
    parfor i=1:size(spec_t,2)
        d=l1*spec_t(:,i);
        lam(:,i)=lsqnonneg(C,d);
    end
    LCs=(vv'*vv)\(vv'*spec_t-D'*lam);
    spec_t_fitted=filter(laser_temp,1,L)*LCs;
    error_temp = spec_t-spec_t_fitted;
    err_shif(si)=norm(error_temp,'fro');
end

[min_err_shif,I_min_err_shif]=min(err_shif);%find best shift
% show result
if show_align_result_flag
    %calcultae best shift
    laser_temp=circshift(laser,[shif(I_min_err_shif),0]);
    vv=filter(laser_temp,1,L);
    cc=vv\spec_t;
    spec_t_fitted=vv*cc;
    decay = L*cc;
    
    h1 = figure('Position',[50 300 1000 400]);
    %     subplot(2,1,1);
    step = ceil(size(spec_t,2)*0.1);
    plot(spec_t(:,1:step:end),'b.--','LineWidth', 1)
    hold on
    %     ylim([-0.05 max(spec_t(:)*1.2)]);
    %     subplot(2,1,2);
    plot(spec_t_fitted(:,1:step:end),'r','LineWidth', 1)
    plot(laser_temp,'g','LineWidth', 1)
    [~, I_max] = max(spec_t(:,1));
    xlim([0 I_max+100])
    ylim([-0.05 max(spec_t(:)*1.2)]);
    legend('Raw Data', 'Fitting','iRF')
    title('iRF alignment result')
    
    h2=figure('Position',[1200 300 500 500]);
    hold on;
    plot(shif,err_shif,'b--','LineWidth',2);
    plot(shif(I_min_err_shif),min_err_shif,'rx','LineWidth',3,'MarkerSize',16);
    % plot(shif(I_min_err_shif-5),err_shif(I_min_err_shif-5),'rx','LineWidth',3,'MarkerSize',16);
    hold off
    pause(0.5)
    
end
out_temp=shif(I_min_err_shif);

% if Running more than 60 sec, then time out
if etime(clock,timeout)> 60
    out = NaN;
    return;
else
%     if I_min_err_shif<=1
%         %     close(h1);
%         out_temp=spec_laser_align(spec,laser,K,n_align_fst,shif_range,shif_c-shif_range/2,timeout);
%         
%     elseif I_min_err_shif>=length(err_shif)-1
%         %     close(h1);
%         out_temp=spec_laser_align(spec,laser,K,n_align_fst,shif_range,shif_c+shif_range/2,timeout);
%     end
    
    out=out_temp;
    if show_align_result_flag
        close(h1)
        close(h2)
    end
end
