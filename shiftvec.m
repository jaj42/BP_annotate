function nvec = shiftvec(vec, offset)
    offset_direction = sign(offset);
    offset_value = abs(offset);
    
    if offset_direction > 0
        nvec = vec(offset_value : end);
        nvec = [nvec zeros(1, offset_value - 1)];
    elseif offset_direction < 0
        nvec = [zeros(1, offset_value) vec];
        nvec = nvec(1 : end - offset_value);
    else
        % offset is 0, no adjustment required.
        nvec = vec;
    end
end