function [L, L_prime, F, F_prime, Bs, Bs_prime, d1, d2] = generate_blocks(N, T, arr_rate, pb, mu_s, mu_b)
   
    % tic
    ps = 1 - pb;
    lambda_s = arr_rate * ps;
    lambda_b = arr_rate * pb;
    d1 = int32(N + ((N * N) / (2 * T)) - (N / 2) + (N / T) + 1);
    d2 = int32((T) * (N / T) + 1);
    [indices_d1, indices_d2] = create_indices(N, T);

    F = zeros(d2, d2);
    for i = 1:size(indices_d2, 1)
        for j = 1:size(indices_d2, 1)
            if isequal(indices_d2(i, :), indices_d2(j, :))
                F(i, j) = arr_rate;
            end
        end
    end

    F_prime = zeros(d1, d2);
    for i = 1:size(indices_d1, 1)
        for j = 1:size(indices_d2, 1)
            if isequal(indices_d1(i, :), indices_d2(j, :))
                if indices_d1(i, 1) + T * (indices_d1(i, 2)) == N
                    F_prime(i, j) = arr_rate;
                end
                if indices_d1(i, 1) + T * (indices_d1(i, 2)) > N - T && indices_d1(i, 1) + T * (indices_d1(i, 2)) < N
                    F_prime(i, j) = lambda_b;
                end
            end
        end
    end

    L_prime = zeros(d1, d1);
    for i = 1:size(indices_d1, 1)
        for j = 1:size(indices_d1, 1)
            if indices_d1(j, 1) == indices_d1(i, 1) + 1 && indices_d1(j, 2) == indices_d1(i, 2) && indices_d1(i, 1) + T * (indices_d1(i, 2)) < N
                L_prime(i, j) = lambda_s;
            end
            if indices_d1(j, 1) == indices_d1(i, 1) && indices_d1(j, 2) == indices_d1(i, 2) + 1 && indices_d1(i, 1) + T * (indices_d1(i, 2)) <= N - T
                L_prime(i, j) = lambda_b;
            end
            if indices_d1(j, 1) == indices_d1(i, 1) - 1 && indices_d1(j, 2) == indices_d1(i, 2)
                L_prime(i, j) = indices_d1(i, 1) * mu_s;
            end
            if indices_d1(j, 1) == indices_d1(i, 1) && indices_d1(j, 2) == indices_d1(i, 2) - 1
                L_prime(i, j) = indices_d1(i, 2) * mu_b;
            end
            if isequal(indices_d1(i, :), indices_d1(j, :))
                L_prime(i, j) = -(1) * (arr_rate + indices_d1(i, 1) * mu_s + indices_d1(i, 2) * mu_b);
            end
        end
    end

    L = zeros(d2, d2);
    for i = 1:size(indices_d2, 1)
        for j = 1:size(indices_d2, 1)
            x = indices_d2(i, :);
            Fr = N - (x(1) + T * x(2));
            if x(1) + T * (x(2)) == N && indices_d2(j, 1) == x(1) - 1 && indices_d2(j, 2) == x(2)
                L(i, j) = x(1) * mu_s * pb;
            end
            if Fr < T - 1 && Fr ~= 0 && indices_d2(j, 1) == x(1) - 1 && indices_d2(j, 2) == x(2)
                L(i, j) = x(1) * mu_s;
            end
            if isequal(indices_d2(i, :), indices_d2(j, :))
                L(i, j) = -(1) * (arr_rate + x(1) * mu_s + x(2) * mu_b);
            end
        end
    end

    Bs_prime = create_Bs_prime_matrices(indices_d2, indices_d1, N, T, mu_s, mu_b, ps, pb, d1, d2);

    Bs = create_B_matrices(indices_d2, N, T, mu_s, mu_b, ps, pb, d1, d2);
    
    % elapsedTime = toc;
    % fprintf('Elapsed time: %.4f Min\n', elapsedTime/60);

    % Find the number of decimal places for each parameter
    arr_rate_decimals = -floor(log10(arr_rate));
    pb_decimals = -floor(log10(pb));
    mu_s_decimals = -floor(log10(mu_s));
    mu_b_decimals = -floor(log10(mu_b));

    % Determine the maximum number of decimal places among the parameters
    max_decimals = max([arr_rate_decimals, pb_decimals, mu_s_decimals, mu_b_decimals]);

    % Construct the format string dynamically
    format_str = ['Blocks-N_%d-T_%d-arrRate_%0.', num2str(max_decimals), 'f-pb_%0.', num2str(max_decimals), 'f-muS_%0.', num2str(max_decimals), 'f-muB_%0.', num2str(max_decimals), 'f.mat'];

    % Use sprintf with the dynamically constructed format string
    mat_name = sprintf(format_str, N, T, arr_rate, pb, mu_s, mu_b);

    %mat_name = sprintf('Blocks-N_%d-T_%d-arrRate_%0.4f-pb_%f-muS_%f-muB_%f.mat', N, T, arr_rate, pb, mu_s, mu_b);
    %save(mat_name, '-v7.3');
end