function [R, pis] = matrix_geom(N, T, arr_rate, pb, mu_s, mu_b, algo, useId, n_pis)

    % Returns matrix R and the vectors of \bold pi_I up to i=n_pis

    % Loading the following package:
    % D. A. Bini, B. Meini, S. Steff ́e, B. Van Houdt, Structured Markov Chains
    % Solver: Algorithms, in: Proceeding from the 2006 Workshop on Tools
    % for Solving Structured Markov Chains, ACM, 2006, p. 13–es.
    
    addpath(fullfile(pwd, '/MG1files'));
    
    
    % Loading functions to construct the matrix blocks
    addpath(fullfile(pwd, '/Aux'));
    
    [L, L_prime, F, F_prime, Bs, Bs_prime, d1, d2] = generate_blocks(N, T, arr_rate, pb, mu_s, mu_b);
    
    l = max([max(max(-L_prime)), max(max(-L))]);
    
    A = [(1/l)*F (1/l)*L+eye(size(L,1))];
    for i=1:length(Bs)
        A = [A (1/l)*Bs{i}];
    end
    
    R = [];
    if useId
        % Using identity matrix;
        [R,numit] = GIM1_R(A,'A', algo, 'StartValue', eye(size(L)));
        disp('Identity matrix');
        disp(numit);
    else
        % 'Using zero matrix;
        [R,numit] = GIM1_R(A,'A', algo);
        disp('Zero matrix');
        disp(numit);
    end
    
    bl = Bs_prime{1};
    for i=2:length(Bs_prime)
        bl = bl + ((R^(i-1))*Bs_prime{i});
    end
    
    br = L;
    for i=1:length(Bs)
        br = br + ((R^i)*Bs{i});
    end
    
    bc = [L_prime F_prime];
    bcb = [bl br];
    bc = [bc; bcb];
    
    bc(1,end) = 1;
    for i=2:size(bc,1)
        bc(i,end) = 0;
    end
    
    a = bc.';
    b = zeros(d1+d2,1);
    b(end) = 1;
    raw_pi = linsolve(a,b);
    
    ir = eye(size(R,1))-R;
    r = raw_pi(d1+1:end).'/ir;
    alfa = sum(raw_pi(1:d1)) + sum(r);
    pi = raw_pi/alfa;
    pis = cell(1,n_pis);
    pis{1} = pi(1:d1).';
    pis{2} = pi(d1+1:end).';
    for i=3:length(pis)
        pis{i} = pis{i-1}*R;
    end

end