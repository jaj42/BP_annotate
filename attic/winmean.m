function [ res ] = winmean( rwin, quant )
    shape = size(rwin);
    vecsize = shape(2);
    res = zeros(1, vecsize);
    for i = 1:vecsize
        res(i) = quant * mean(rwin(:,i));
    end
end