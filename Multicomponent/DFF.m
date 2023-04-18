function [val] = DFF(phi,N)
%derivative of F is, which is a vector
    val = zeros(size(phi));
    for i=1:N
        for j=1:N
            if i~=j
                val(:,i) = val(:,i) + 4*phi(:,i).*phi(:,j).^2;
            end
        end
    end   
end


