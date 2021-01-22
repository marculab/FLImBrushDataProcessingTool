function Background = LoadBackground(filepath)
    fid = fopen(filepath,'r');
    if fid < 0
        error('Cannot open background file "%s"!',filepath);
    else
        head = fread(fid,[1 4],'int32');
        Background = fread(fid,[1 head(4)],'double');
        fclose(fid);
    end

end