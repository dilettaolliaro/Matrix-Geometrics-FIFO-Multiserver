% Here we perform an explanotory experiment in particular we compute
% performance metrics for different values of big and small probabilties

% In order to be able to compare these results we will fix the desired load
% factor rho for the system and compute the arrival rate lambda (l) for 
% each specific configuration accordingly
syms l

% Number of Available Servers in the Systems
N = 200;
% Amount of Servers Required by a Job Belonging to the Big Class
T = 100;
% Small Class Jobs Service Time
taus = 1;
% Small Class Jobs Service Rate
mu_s = 1/taus;
% Big Class Jobs Service Time
taub = 4;
% Big Class Jobs Service Time
mu_b = 1/taub;
% Probability of Small Class
ps = 0.8;
% Probability of Big Class
pb = 0.2;

% Values of rho we are going to compute metrics for
rhos = linspace(0.5, 0.8, 30);
 
filename = sprintf('Results/experiment_overLoad_N%d_T%d_ps%.2f_mu_s%.2f_mu_b%.2f.csv', N, T, ps, mu_s, mu_b);

for i=1:numel(rhos)

    rho = rhos(i);
    
    eq = rho - ((l*(ps*taus + pb * taub * T))/N) == 0;
    arr_rate = double(vpasolve(eq, l)); 
            
    [Nt,Ns,Nw,Ntt,Nst,Nwt,RT,RTs,RTb,WT,WTs,WTb,U,B,Thr,Thrs,Thrb,Wasted,Wasted_HOL,Wasted_no_queue,exeTime] = performanceMetrics(N, T, arr_rate, pb, mu_s, mu_b);   
        
    % Open or create the CSV file in append mode
    fid = fopen(filename, 'a');

    % Write headers if the file is empty
    if ftell(fid) == 0
        fprintf(fid, ['Load,Arr. Rate,T1 Queue,T%d Queue,Queue Total,T1 Service,T%d Service,'...
            'T1 System,T%d System,Avg. System,T1 Waiting,T%d Waiting,WaitTime Total,' ...
            'T1 RespTime,T%d RespTime,RespTime Total,T1 Throughput,T%d Throughput,Utilization,Avg. Busy,' ...
            'Wasted-Tot,Wasted-HoL,Wasted-EmptyQ,Exe. Time\n'], [T, T, T, T, T, T]);
    end

    % Write data to the file
    fprintf(fid, '%f, %f, %f,  %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %.4f\n',...
        rho, arr_rate, Nw(1), Nw(2), Nwt, Ns(1), Ns(2), Nt(1), Nt(2), Ntt, WTs, WTb, WT, RTs, RTb, RT, Thrs, Thrb, U, B, Wasted, Wasted_HOL, Wasted_no_queue, exeTime);

    % Close the file
    fclose(fid);
end












