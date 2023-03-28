function [L2] = calcNorm(err,M,number_of_timesteps)
    L2 = zeros(1,number_of_timesteps);
    for i=1:number_of_timesteps
        L2(i) = sqrt(transpose(err(:,i))*M*err(:,i));
    end
end

