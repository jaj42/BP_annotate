function [ rwin ] = rollingWindow( vector, winsize )
    vector = vector(:);
    vecsize = length(vector);
    rwin = NaN(winsize, vecsize);
    for i = 1 : winsize
        tmp = vector(1 : end - i + 1);
        rwin(i, i : end) = tmp;
    end
end