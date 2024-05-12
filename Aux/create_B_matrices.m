function Bs = create_B_matrices(indices_d2, N, T, mu_s, mu_b, ps, pb, d1, d2)
    Bs = cell(1, T);
    B_1 = zeros(d2, d2);
    for i = 1:size(indices_d2, 1)
        x = indices_d2(i, :);
        for j = 1:size(indices_d2, 1)
            y = indices_d2(j, :);
            if x(1) + T * x(2) == N
                if isequal(x, y)
                    B_1(i, j) = x(1) * mu_s * ps;
                elseif y(1) == x(1) + 1 && y(2) == x(2) - 1
                    B_1(i, j) = x(2) * mu_b * pb * ps;
                end
            end
            Fr = N - (x(1) + T * x(2));
            if Fr == T - 1
                if y(1) == x(1) - 1 && y(2) == x(2) + 1
                    B_1(i, j) = x(1) * mu_s;
                end
            end
            if isequal(x, y)
                B_1(i, j) = B_1(i, j) + x(2) * mu_b * pb;
            end
        end
    end
    Bs{1} = B_1;

    for l = 2:T
        B_l = zeros(d2, d2);
        for i = 1:size(indices_d2, 1)
            x = indices_d2(i, :);
            for j = 1:size(indices_d2, 1)
                y = indices_d2(j, :);
                if x(1) + T * x(2) == N
                    if y(1) == x(1) + l && y(2) == x(2) - 1
                        B_l(i, j) = x(2) * mu_b * (ps ^ l) * pb;
                    end
                    if y(1) == x(1) + T && y(2) == x(2) - 1 && T == l
                        B_l(i, j) = x(2) * mu_b * (ps ^ T);
                    end
                end
                Fr = N - (x(1) + T * x(2));
                if Fr <= T - 1 && Fr ~= 0
                    if y(1) == x(1) + (l - 1) && y(2) == x(2)
                        B_l(i, j) = x(2) * mu_b * (ps ^ (l - 1)) * pb;
                    end
                    if y(1) == x(1) + Fr && y(2) == x(2) && l - 1 == Fr
                        B_l(i, j) = x(2) * mu_b * (ps ^ Fr);
                    end
                end
            end
        end
        Bs{l} = B_l;
    end
end