function [A,T,avglife,intensity,fitt,raw]=monoexp_fit(spec,dt,laser)
A =[];T=[];h =[];fitt=[];raw=[];
%figure
parfor ii = 1:size(spec,2)
    if any((spec(:,ii))>0)
        Y = spec(:,ii)./max(spec(:,ii));
        lower = [0 0.01];
        upper = [1 40];
        %lower = [0 0 eps eps];
        %upper = [1 1 Inf Inf];
        [~,b] = max(Y);
        ynew = Y(b:end);
        if length(ynew)>3
            tt = linspace(0,length(ynew)-1,length(ynew))*dt;
%             fitWeights = ynew>0.2; %set weights of the fitting to 0 for tails
            op = fitoptions('Method', 'NonlinearLeastSquares','Lower',lower,'Upper',upper);
%             op.Weights = fitWeights;
            ft2 = fittype('monoexp_model_init(x,a,t)','options',op);
            %ft2 = fittype('biexp_model_init2(x,a1,a2,t1,t2)','options',op);
            
            % suppress warning of no start point
            warning('off','curvefit:fit:noStartPoint')
            [ff,~,~] = fit(tt',ynew,ft2);
            
            %     subplot(1,2,1),plot(ynew)
            %     decay = ff.a1.*exp(-tt./ff.t1)+(1-ff.a1).*exp(-tt./ff.t2);
            %     hold on
            %     plot(decay,'r'), legend('raw','fit'), title('initial fit')
            %     hold off
            start = [ff.a ff.t];
        else
            start = [0.5 1];
        end
        %start = [ff.a1 ff.a2 ff.t1 ff.t2];
        x = linspace(0,length(Y)-1,length(Y))*dt; % time in ns
        op = fitoptions('Method', 'NonlinearLeastSquares','Lower',lower,'Upper',upper,'StartPoint',start);
        fitWeights2 = zeros(size(Y));
        fitWeights2(1:90) = 1;
        op.Weights = fitWeights2;
        ft = fittype('monoexp_model(x,a,t,L)','problem','L','options',op);
        %ft = fittype('biexp_model3(x,a1,a2,t1,t2,L)','problem','L','options',op);
        
        [f,gof,gg] = fit(x',Y,ft,'problem',laser);
        
        decay = f.a.*exp(-x./f.t);
        %decay = f.a1.*exp(-x./f.t1)+(f.a2).*exp(-x./f.t2);
%         factor = sum((spec(:,ii)))/sum(decay);
%         decay = decay*factor;
        y = filter(laser,1,decay);% fitting
               
        decay = decay/max(y)*max(spec(:,ii));
        yy = filter(laser,1,decay); 
        
        y = y./max(y); 
        fff = y(1:length(x));
        yyy = Y;
        %plot normalized fitting result
%         plot(x, yy, 'LineWidth', 2)
%         hold on
%         plot(x, spec(:,ii), 'r-'),legend('fit','raw'), title([num2str(ii) ' and ' num2str(f.t)]),hold off
%         pause(.1)
        
        A = [A f.a];
        %A2 = [A2 f.a2];
        T=[T f.t];
        h = [h,decay'];
        fitt = [fitt,fff'];
        raw = [raw,yyy];
    else
        A = [A 0];
        T = [T NaN];
        h = [h, spec(:,ii)];
        fitt = [fitt,spec(:,ii)];
        raw = [raw,spec(:,ii)];
    end
end
%     Ts = [T1;T2];As=[A1;A2];
%     [TF, I] = sort(Ts,1);
%     AF = zeros(size(As));
%     for k =1: size(I,2)
%         AF(:,k) = As(I(:,k),k);
%     end
%  T1 = TF(1,:);T2=TF(2,:);A1=AF(1,:);A2=AF(2,:);
% average lifetime from decay
[avglife,intensity]=h_lifet(h,dt,'average');

% re-enable the warning
warning('on','curvefit:fit:noStartPoint')
end