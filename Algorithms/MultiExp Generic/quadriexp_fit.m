function [A1,A2,A3,A4,T1,T2,T3,T4,avglife,intensity,fitt,raw,h]=quadriexp_fit(spec,dt,laser)
A1 =[];A2=[];A3=[];A4=[];T1=[];T2=[];T3=[];T4=[];h =[];fitt=[];raw=[];
%figure
parfor ii = 1:size(spec,2)
    if any((spec(:,ii))>0)
        Y = spec(:,ii);
        lower = [0 0 0 0 0.01 0.01 0.01 0.01];
        upper = [1 1 1 1 25 25 25 25];
        %lower = [0 0 eps eps];
        %upper = [1 1 Inf Inf];
        [~,b] = max(Y);
        ynew = Y(b:end);
        %     if ~any(isnan(ynew(:)))
        if length(ynew)>3
            tt = linspace(0,length(ynew)-1,length(ynew))*dt;
            op = fitoptions('Method', 'NonlinearLeastSquares','Lower',lower,'Upper',upper);
            ft2 = fittype('quadriexp_model_init(x,a1,a2,a3,a4,t1,t2,t3,t4)','options',op);

            % suppress warning of no start point
            warning('off','curvefit:fit:noStartPoint')
            [ff,~,~] = fit(tt',ynew,ft2);
            
            %     subplot(1,2,1),plot(ynew)
            %     decay = ff.a1.*exp(-tt./ff.t1)+(1-ff.a1).*exp(-tt./ff.t2);
            %     hold on
            %     plot(decay,'r'), legend('raw','fit'), title('initial fit')
            %     hold off
            start = [ff.a1 ff.a2 ff.a3 ff.a4 ff.t1 ff.t2 ff.t3 ff.t4];
        else
            start = [0.5 0.5 0.5 0.5 1 1 1 1];
        end
        %start = [ff.a1 ff.a2 ff.t1 ff.t2];
        x = linspace(0,length(Y)-1,length(Y))*dt; % time in ns
        op = fitoptions('Method', 'NonlinearLeastSquares','Lower',lower,'Upper',upper,'StartPoint',start);
        ft = fittype('quadriexp_model(x,a1,a2,a3,a4,t1,t2,t3,t4,L)','problem','L','options',op);
        [f,gof,gg] = fit(x',Y,ft,'problem',laser);
        decay = f.a1.*exp(-x./f.t1)+f.a2.*exp(-x./f.t2)+f.a3.*exp(-x./f.t3)+f.a4.*exp(-x./f.t4);
        y = filter(laser,1,decay);
        fff = y;
        yyy = Y;
        %     plot(x, y(1:length(x)), 'LineWidth', 2)
        %     hold on
        %     plot(x, Y, 'r-'),legend('fit','raw'), title(num2str(ii)),hold off
        %     pause(.1)
        A1 = [A1 f.a1];
        A2 = [A2 f.a2];
        A3 = [A3 f.a3];
        A4 = [A4 f.a4];
        T1 = [T1 f.t1];
        T2 = [T2 f.t2];
        T3 = [T3 f.t3];
        T4 = [T4 f.t4];
        h = [h,decay'];
        fitt = [fitt,fff'];
        raw = [raw,yyy];
    else
        A1 = [A1 0];
        A2 = [A2 0];
        A3 = [A3 0];
        A4 = [A4 1];
        T1 = [T1 NaN];
        T2 = [T2 NaN];
        T3 = [T3 NaN];
        T4 = [T4 NaN];
        h = [h,zeros(size(spec(:,ii)))];
        fitt = [fitt,spec(:,ii)];
        raw = [raw,spec(:,ii)];
    end
end
Ts = [T1;T2;T3;T4];As=[A1;A2;A3;A4];
[TF, I] = sort(Ts,1);
AF = zeros(size(As));
for k =1: size(I,2)
    AF(:,k) = As(I(:,k),k);
end
T1 = TF(1,:);T2=TF(2,:);T3 = TF(3,:);T4 = TF(4,:);A1=AF(1,:);A2=AF(2,:);A3=AF(3,:);A4=AF(4,:);
% average lifetime from decay
[avglife,intensity]=h_lifet(h,dt,'average');

% re-enable the warning
warning('on','curvefit:fit:noStartPoint')
end