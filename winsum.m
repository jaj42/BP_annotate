function [ res ] = winsum( rwin )
%WINSUM Summary of this function goes here
%   Detailed explanation goes here
    shape = size(rwin);
    if length(shape) > 2
        error('Incorrect window');
    end
    
    vecsize = shape(2);
    res = zeros(1, vecsize);
    
    for i = 1:vecsize
        res(:,i) = sum(rwin(:,i));
    end
end