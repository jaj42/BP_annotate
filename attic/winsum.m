function [ res ] = winsum( rwin )

    shape = size(rwin);
    vecsize = shape(2);
    res = zeros(1, vecsize);
    for i = 1:vecsize
        res(i) = sum(rwin(:,i));
    end
end