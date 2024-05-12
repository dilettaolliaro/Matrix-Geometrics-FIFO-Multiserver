function Bs_prime = create_Bs_prime_matrices(indices_d2, indices_d1, N, T, mu_s, mu_b, ps, pb, d1, d2)
    Bs_prime = cell(1, T);
    B1_prime = zeros(d2, d1);
    for i = 1:size(indices_d2, 1)
        x = indices_d2(i, :);
        Fr = N - (x(1) + T * x(2));
        for j = 1:size(indices_d1, 1)
            y = indices_d1(j, :);
            if x(1) + T * x(2) == N
                if isequal(x, y)
                    B1_prime(i, j) = x(1) * mu_s * ps + x(2) * mu_b * pb;
                elseif y(1) == x(1) + 1 && y(2) == x(2) - 1
                    B1_prime(i, j) = x(2) * mu_b * ps;
                end
            end
            if Fr ~= 0 && Fr <= T - 1
                if isequal(x, y)
                    B1_prime(i, j) = x(2) * mu_b;
                end
            end
            if Fr == T - 1
                if y(1) == x(1) - 1 && y(2) == x(2) + 1
                    B1_prime(i, j) = x(1) * mu_s;
                end
            end
        end
    end
    Bs_prime{1} = B1_prime;

    for l = 2:T
        Bl_prime = zeros(d2, d1);
        for i = 1:size(indices_d2, 1)
            x = indices_d2(i, :);
            Fr = N - (x(1) + T * x(2));
            for j = 1:size(indices_d1, 1)
                y = indices_d1(j, :);
                if x(1) + T * x(2) == N && y(1) == x(1) + l && y(2) == x(2) - 1
                    Bl_prime(i, j) = x(2) * mu_b * (ps ^ l);
                end
                if Fr <= T - 1 && Fr ~= 0 && y(1) == x(1) + (l - 1) && y(2) == x(2)
                    Bl_prime(i, j) = x(2) * mu_b * (ps ^ (l - 1));
                end
            end
        end
        Bs_prime{l} = Bl_prime;
    end
end