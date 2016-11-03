function [ rwin ] = rollingWindow( vector, winsize )
%ROLLINGWINDOW Return a rolling window with aperture set to winsize
    vector = vector(:);
    vecsize = length(vector);
    buffer = NaN(winsize, vecsize);

    for i = 1 : winsize
        tmp = vector(1 : end - i + 1);
        buffer(i, i : end) = tmp;
    end
    
    rwin = buffer;
end