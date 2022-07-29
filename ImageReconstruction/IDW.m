%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INVERSE DISTANCE WEIGHT %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Vint]=IDW(xc,yc,vc,x,y,e,r1,r2,SNR)
%%%%%%%%%%%%%%%%%
%%% INPUTS
%xc = stations x coordinates (columns) [vector]
%yc = stations y coordinates (rows) [vector]
%vc = variable values on the point [xc yc]
%x = interpolation points  x coordinates [vector]
%y = interpolation points y coordinates [vector]
%e = distance weight
%r1 --- 'fr' = fixed radius ;  'ng' = neighbours
%r2 --- radius lenght if r1 == 'fr' / number of neighbours if  r1 =='ng'
%%% OUTPUTS
%Vint --- Matrix [length(y),length(x)] with interpolated  variable values
%%% EXAMPLES
%%% --> V_spa=IDW(x1,y1,v1,x,y,-2,'ng',length(x1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vint=zeros(length(x),1);
% xc=reshape(xc,1,length(xc));
% yc=reshape(yc,1,length(yc));
% vc=reshape(vc,1,length(vc));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  strcmp(r1,'fr')
    if  (r2<=0)
        disp('Error: Radius must be positive')
        return
    end
    for i=1:length(x)
%         if (x(i)==958)&&(y(i)==582)
%             disp('Pause')
%         end
        
        D=[]; V=[]; wD =[]; vcc=[]; wSNR = []; w = [];
        D= sqrt((x(i)-xc).^2 +(y(i)-yc).^2);
        if min(D)==0 % point is a measurement point
            IdxX = find(x(i)==xc);
            IdxY = find(y(i)==yc);
            Idx = intersect(IdxX,IdxY);
            [~,I] = max(SNR(Idx)); %find the point with max SNR
            Vint(i)=vc(Idx(I)); % if there are duplicate, use the 1st one
        else
            vcc=vc(D<r2); SNR_temp = SNR(D<r2); D=D(D<r2);
%             wD = ((r2-D)/r2).^e;
            wD = D.^e;
            wSNR = SNR_temp./sum(SNR_temp);
            w = wD.*wSNR;
            V = vcc.*w;
            
            if isempty(D)
                V=NaN;
            else
                V=sum(V)/sum(w);
            end
            Vint(i)=V;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    if (r2 > length(vc)) || (r2<1)
        disp('Error: Number of neighbours not congruent with data')
        return
    end
    for i=1:length(x)
        D=[]; V=[]; wD =[];vcc=[];
        D= sqrt((x(i)-xc).^2 +(y(i)-yc).^2);
        if min(D)==0
            disp('Error: One or more stations have the coordinates of an interpolation point')
            return
        end
        [D,I]=sort(D);
        vcc=vc(I);
        V = vcc(1:r2).*(D(1:r2).^e);
        wD = D(1:r2).^e;
        V=sum(V)/sum(wD);
        Vint(i)=V;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
