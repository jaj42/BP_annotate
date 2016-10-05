function [ res ] = winquant( rwin, quant, varargin )
%WINQUANT Summary of this function goes here
%   Detailed explanation goes here
    shape = size(rwin);
    if length(shape) > 2
        error('Incorrect window');
    end
    
    if nargin < 2
        quant = .7;
    end
    
    vecsize = shape(2);
    res = zeros(1, vecsize);

    for i = 1:vecsize
        res(:,i) = quantile(rwin(:,i), quant);
    end
end