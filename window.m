function [ rwin ] = window( vector, winsize )
%WINDOW Summary of this function goes here
%   Detailed explanation goes here
    vector = vector(:);
    vecsize = length(vector);
    %buffer = zeros(winsize, vecsize);
    buffer = NaN(winsize, vecsize);
    buffer(1, :) = vector;
    for i = 1 : winsize
        tmp = vector(i : end);
        buffer(i, 1 : length(tmp)) = tmp;
    end
    rwin = buffer;
end
