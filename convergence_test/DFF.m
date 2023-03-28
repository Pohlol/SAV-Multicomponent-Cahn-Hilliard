function [val] = DFF(phi,N)
    val = zeros(size(phi));
    for i=1:N
        for j=1:N
            if i~=j
                val(:,i) = val(:,i) + 2*phi(:,i).*phi(:,j).^2;
            end
        end
        val(:,i) = phi(:,i);
    end   

end


