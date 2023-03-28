function [val] = projecton(DF,N)
    %projection on the space where the vector elements add up to 0
    val = DF - 1/N*sum(DF,3);
end

