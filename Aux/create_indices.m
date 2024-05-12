function [indices_d1, indices_d2] = create_indices(N, T)
    indices_d1 = [];
    indices_d2 = [];
    for s = 0:N
        for b = 0:round(N/T)
            if s + (b * T) <= N
                indices_d1 = [indices_d1; [s, b]];
            end
            if s + (b * T) > (N - T) && s + (b * T) <= N
                indices_d2 = [indices_d2; [s, b]];
            end
        end
    end
end