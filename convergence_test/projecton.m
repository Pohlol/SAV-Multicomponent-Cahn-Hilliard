function [val] = projecton(DF,N)
    val = DF - 1/N*sum(DF,3);
end

