function exponentials = multiexp_model(A, T, timebase, varargin)

% A - Weights
% T - Taus (in ns)
% timebase - timebase of the exponential functions (in ns)
% intensity - overall intensity
timebase = reshape(timebase, length(timebase),1);
orders = size(T,1);
points = size(T,2);
exponentials = cell(points,1);


% if intensity is given
if nargin > 3
    intensity = varargin{1};
    scaled = 1;
else
    intensity = ones(points,1);
    scaled = 0;
end

for i = 1:points
    expfunc = zeros(length(timebase),orders);
    for j = 1:orders
        expfunc(:,j) = A(j,i).*exp(-timebase./T(j,i));
    end
%     expfunc=expfunc';
    if scaled
        exponentials{i} = intensity(i)* expfunc ./ sum(expfunc(:));
    else
        exponentials{i} = expfunc ./ sum(expfunc(:));
% diable scaling
%         exponentials{i} = expfunc;
    end
    
end
end