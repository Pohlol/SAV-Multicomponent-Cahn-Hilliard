function [val] = FF(phi,N)
    %polynomial Potential for Multicomponent Cahn Hilliard Equations with a
    %small positive constant added
    val = 0.0001;
    for i=1:N
        for j=1:N
            if i~=j
                val = val + phi(:,i).^2.*phi(:,j).^2;
            end
        end
    end
end

