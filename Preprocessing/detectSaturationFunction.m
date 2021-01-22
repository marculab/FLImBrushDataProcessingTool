function dataInd = detectSaturation(spec, minVal, maxVal)

data = spec;

dataInd = zeros(size(data, 2), 1);
k = 1;
%loop through all data points
for i = 1:size(data, 2)
    
    if ge(max(data(:, i)), maxVal) % check saturation
        
        warning('off')
        [pks, ~] = find(ge(data(:, i), (maxVal - eps)));
        
        if and(~isempty(pks), ge(size(pks, 1), 2))
            
            if lt(abs(min(diff(pks))), 3)
                
                dataInd(k) = i;
                k = k + 1;
                
            end
         
        end
        
    elseif le(max(data(:, i)), minVal) % check low data
        
        dataInd(k) = i;
        k = k + 1;
        
    end

end

dataInd(k:end) = []; % delect empty element