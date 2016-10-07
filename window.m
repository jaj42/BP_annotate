function [ rwin ] = window( vector, winsize )
%WINDOW Summary of this function goes here
%   Detailed explanation goes here
    vector = vector(:);
    vecsize = length(vector);
    buffer = NaN(winsize, vecsize);

    for i = 1 : winsize
        tmp = vector(1 : end - i + 1);
        buffer(i, i : end) = tmp;
    end
    
    rwin = buffer;
end