function [A1,A2,T1,T2,avglife,intensity,fitt,raw,h]=biexp_fit(spec,dt,laser,taus)
A1 =[];A2=[];T1=[];T2=[];h =[];fitt=[];raw=[];

if isempty(taus) % if taus is empty, use defaut
    lower = [0 0 0.01 0.01];
    upper = [Inf Inf 25 25];
else % if taus is not empty, fix tau
    lower = [0 0 taus(1) taus(2)];
    upper = [Inf Inf taus(1) taus(2)];
end

parfor ii = 1:size(spec,2)
    if any((spec(:,ii))>0)
        Y = spec(:,ii);% ./max(spec(:,ii));
        
        [~,b] = max(Y);
        ynew = Y(b:end);
        %         Y=ynew;
        if length(ynew)>3
            tt = linspace(0,length(ynew)-1,length(ynew))*dt;
            op = fitoptions('Method', 'NonlinearLeastSquares','Lower',lower,'Upper',upper);
            ft2 = fittype('biexp_model_init(x,a1,a2,t1,t2)','options',op);
            %ft2 = fittype('biexp_model_init2(x,a1,a2,t1,t2)','options',op);
            
            % suppress warning of no start point
            warning('off','curvefit:fit:noStartPoint')
            [ff,~,~] = fit(tt',ynew,ft2);
            
            %     subplot(1,2,1),plot(ynew)
            %     decay = ff.a1.*exp(-tt./ff.t1)+(1-ff.a1).*exp(-tt./ff.t2);
            %     hold on
            %     plot(decay,'r'), legend('raw','fit'), title('initial fit')
            %     hold off
            start = [ff.a1 ff.a2 ff.t1 ff.t2];
            %             start = [0.5 1 1];
        else
            start = [0.5 0.5 2 2]; % starting points
        end
        x = linspace(0,length(Y)-1,length(Y))*dt; % time in ns
        op = fitoptions('Method', 'NonlinearLeastSquares','Lower',lower,'Upper',upper,'StartPoint',start);
        ft = fittype('biexp_model(x,a1,a2,t1,t2,L)','problem','L','options',op);
        %ft = fittype('biexp_model3(x,a1,a2,t1,t2,L)','problem','L','options',op);
        [f,gof,gg] = fit(x',Y,ft,'problem',laser);
        
        %         plot(x',Y);
        %         hold on
        decay = f.a1.*exp(-x./f.t1)+f.a2.*exp(-x./f.t2);
        %decay = f.a1.*exp(-x./f.t1)+(f.a2).*exp(-x./f.t2);
        % correct for intensity
        %         factor = sum((spec(:,ii)))/sum(decay);
        %         decay = decay*factor;
        y = filter(laser,1,decay);
        %         y = y./max(y);
        fff = y;
        yyy = Y;
        %         plot(x,fff);
        %         legend('off');
        
        %     plot(x, y(1:length(x)), 'LineWidth', 2)
        %     hold on
        %     plot(x, Y, 'r-'),legend('fit','raw'), title(num2str(ii)),hold off
        %     pause(.1)
        
        A1 = [A1 f.a1];T1=[T1 f.t1];T2=[T2 f.t2];
        A2 = [A2 f.a2];
        %A2 = [A2 f.a2];
        h = [h,decay'];
        fitt = [fitt,fff'];
        raw = [raw,yyy];
    else
        A1 = [A1 0];
        A2 = [A2 1];
        T1=[T1 NaN];
        T2=[T2 NaN];
        h = [h, spec(:,ii)];
        fitt = [fitt,spec(:,ii)];
        raw = [raw,spec(:,ii)];
    end
    
end
% close(hh)
% reorder according to lifetime values
Ts = [T1;T2];As=[A1;A2];
[TF, I] = sort(Ts,1);
AF = zeros(size(As));
for k =1: size(I,2)
    AF(:,k) = As(I(:,k),k);
end
T1 = TF(1,:);T2=TF(2,:);A1=AF(1,:);A2=AF(2,:);
% A1 = A1.*MaxData;
% A2 = A2.*MaxData;
% average lifetime from decay
intensity = sum(h);
avglife = (A1.*T1.^2+A2.*T2.^2)./(A1.*T1+A2.*T2);
% [avglife,intensity]=h_lifet(h,dt,'average');
% intensity = Raw_INT;
% re-enable the warning
warning('on','curvefit:fit:noStartPoint')
end