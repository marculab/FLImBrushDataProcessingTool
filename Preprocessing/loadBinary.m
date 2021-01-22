function [varargout]=loadBinary(filepath,type,saveflag)
% LOADBINARY load binary data saved by Labview
%   data = loadBinary(filepath,type,saveflag)
%   filepath - string of the binary data path/filename.
%   type - string of data type. e.g. 'int8', 'int16', 'double'
%   saveflag - bool flag to choose saving/not saving loaded data in .mat
%       format with the SAME name
%   Create and updated by Michael Ma by June 2014
%   
fid = fopen(filepath);
if fid ~= -1
    n = fread(fid,[1,1],'int32');
    m = fread(fid,[1,1],'int32');
    idx=1;
    N = 0;
    while ~feof(fid)
        newdata = fread(fid,[m,n],type);
        data2{idx} = newdata';
        idx = idx + 1;
        N = N + n;
        n = fread(fid,[1,1],'int32');
        m = fread(fid,[1,1],'int32');
    end
    chunks = length(data2);
    m=size(data2{1},2);
    n=size(data2{1},1);
    alldata = zeros(N,m);
    firstlines = zeros(chunks,m);
    fileptr = 0;
    for idx = 1:chunks
        chunksize=size(data2{idx},1);
        alldata((1:chunksize)+fileptr,:) = data2{idx};
        firstlines(idx,:) = data2{idx}(1,:);
        fileptr = fileptr + chunksize;
    end
    if ~~saveflag
        save([filepath '.mat'],'data');
    end
    fclose(fid);
    
    varargout{1} = alldata;
    if nargout > 1    
        varargout{2} = firstlines;
    end
else
    error('The file "%s" cannot be opened!',filepath);
end

end