function gain = get_ula_gain(N,w,az)
    array_response = get_ula_response(N, az);

    gain = conj(w)' * array_response;
end

