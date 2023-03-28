function [val] = FF(phi,N)
    val = 0.0001;
    for i=1:N
        for j=1:N
            if i~=j
                val = val + phi(:,i).^2.*phi(:,j).^2;
            end
        end
        val = val + 1/2*phi(:,i).^2;
    end
end

