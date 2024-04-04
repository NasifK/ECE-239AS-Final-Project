function [response_vector] = get_ula_response(N,az)
    % assume az is in degrees 
    az = deg2rad(az);
    response_vector = zeros(N,1);
    for i = 1:N
        response_vector(i) = exp(-1i * pi * i * sin(az));
    end
end

