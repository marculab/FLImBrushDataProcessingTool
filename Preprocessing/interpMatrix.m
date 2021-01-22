function interpedData = interpMatrix(Sample,data,query)
interpedData = data;
for i = 1:size(Sample,1)
    interpedData(i,:) = interp1(Sample(i,:),data(i,:),query);
end