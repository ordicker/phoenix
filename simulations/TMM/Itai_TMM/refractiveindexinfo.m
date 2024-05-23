
function [epss] = refractiveindexinfo(filename,freq)
    c=299792458; %m/s
    lines = regexp(fileread(filename),'\n','split');

    data=readmatrix(filename);
    whichline = find(isnan(data(:,1))); % find wl,k line
    wl_n = data(1:whichline-1,:);
    wl_k = data(whichline+1:end,:);
    
    n = interp1(wl_n(:,1), wl_n(:,2),c./freq*1e6);
    k = interp1(wl_k(:,1), wl_k(:,2),c./freq*1e6);
    epss = (n+1i*k).^2
    
end
