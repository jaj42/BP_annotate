function [ res ] = winquant( rwin, quant )
%WINQUANT Summary of this function goes here
%   Detailed explanation goes here
    shape = size(rwin);
    vecsize = shape(2);

    res = zeros(1, vecsize);

    for i = 1:vecsize
        res(i) = quantile(rwin(:,i), quant);
    end
end