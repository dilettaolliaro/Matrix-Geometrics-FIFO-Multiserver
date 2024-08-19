function [Nt,Ns,Nw,Ntt,Nst,Nwt,RT,RTs,RTb,WT,WTs,WTb,U,B,Thr,Thrs,Thrb,Wasted,Wasted_busy_queue,Wasted_no_queue,exeTime] = performanceMetrics(N, T, arr_rate, pb, mu_s, mu_b)

tic
ps = 1-pb;

% change algo to one of the following if needed:
%       'CR' : Cyclic Reduction [Bini, Meini]
%       'NI' : Newton Iteration [Perez, Telek, Van Houdt]
%       'RR' : Ramaswami Reduction [Bini, Meini, Ramaswami]
%       'IS' : Invariant Subspace [Akar, Sohraby]

% Default is: 'FI' : Functional Iterations [Neuts]
algo = 'FI';
% set to 1 to start R computation from 0 matrix or from the identity matrix
% By default we start fromt he identity matrix
useId = 1;
util = (arr_rate*((pb*(1/mu_b)*T)+(ps*(1/mu_s)*1)))/N;
if util <= 0.3
    useId = 0;
end
% n_pis gives the number i of \bold pi_i vectors produced by matrix_geom
% function
n_pis = 100;

% Returning:
    % Nt = [Nt_s, Nt_b] = # [small, big] jobs in system
    % Ns = [Ns_s, Ns_b] = # [small, big] jobs in service
    % Nw = [Nw_s, Nw_b] = # [small, big] jobs in waiting queue
    % Ntt = total # jobs in system
    % Nst = total # jobs in service
    % Nwt = total # jobs in waiting queue
    % RT = Overall Avg. Resp. Time
    % RTs = Small Class Avg. Resp. Time
    % RTb = Big Class Avg. Resp. Time
    % WT = Overall Avg. Wait. Time
    % WTs = Small Class Avg. Wait. Time
    % WTb = Big Class Avg. Wait. Time
    % U = Utilization
    % B = Avg. Number of Busy Servers
    % Thr = Overall System Throughput
    % Thrs = Small Class Throughput
    % Thrb = Big Class Throughput
    % Wasted = Avg. Number of Idle Servers
    % Wasted_busy_queue = Avg. Number of Wasted Servers (because of HOL)
    % Wasted_no_queue = Avg. Number of Idle Servers (because empty Queue)
    % exeTime = Execution Time


Lz = [];
for s=0:N
    for b=0:N
        if s+b*T <= N
            Lz = [Lz; s b];
        end
    end
end

L = [];
H = [];
Lw = [];
for s=0:N
    for b=0:N
        if s+b*T <= N && s+b*T > (N-T)
            if s+b*T < N
                L = [L; s b+1];
                Lw = [Lw; 0 1];
            else
                L = [L; s+ps b+pb];
                Lw = [Lw; ps pb];
            end
            H = [H; s b];
        end
    end
end

% Compuute matrix R and \bold pi_i vectors up to i = n_pis

[R, pis] = matrix_geom(N, T, arr_rate, pb, mu_s, mu_b, algo, useId, n_pis);

% N
Nl = pis{1}*Lz;
Nm = [ps pb]*sum((pis{2}*R*((eye(size(R,1))-R)^-2)));
Nr = pis{2}*(inv(eye(size(R,1))-R))*L;

% Average Number of Jobs in the System - Vector of two [small and big]
Nt = Nl+Nm+Nr;

% Average Total Number of Jobs in the System
Ntt = sum(Nt);

% Average Number of Jobs Service - Vector of two [small and big]
Ns = pis{1}*Lz + pis{2}*(inv(eye(size(R,1))-R))*H;

% Average Total Number of Jobs in Service
Nst = sum(Ns);

% Average Number of Jobs in the Waiting Line - Vector of two [small and big]
Nw = [ps pb]*sum((pis{2}*R*((eye(size(R,1))-R)^-2))) + pis{2}*(inv(eye(size(R,1))-R))*Lw;

% Average Total Number of Jobs in the Waiting Line
Nwt = sum(Nw);

% Average Number of Busy Servers
B = pis{1}*Lz*[1 T].' + pis{2}*inv(eye(size(R,1))-R)*H*[1 T].';
% Utilisation (i.e. Average Number of Busy Servers over Available Servers)
U = B/N;


% Average Response Time for the Overall System
RT = sum(Nt)/arr_rate;
% Average Response Time for the Small Job Class 
RTs = Nt(1)/(arr_rate*ps);
% Average Response Time for the Big Job Class 
RTb = Nt(2)/(arr_rate*pb);


% Average Wating Time for the Overall System
WT = sum(Nw)/arr_rate;
% Average Wating Time for the Small Job Class 
WTs = Nw(1)/(arr_rate*ps);
% Average Wating Time for the Big Job Class 
WTb = Nw(2)/(arr_rate*pb);


% Overall System Throughput
Thr = sum(Nt)/RT;
% Small Class Throughput
Thrs = Nt(1)/RTs;
% Big Class Throughput
Thrb = Nt(2)/RTb;


% We disaggregate the number of idle servers, distinguishing from those idle 
% because of an empty queue and those idle EVEN IF the queue is not empty 
% (the latter are going to be addressed as Wasted)

% Probability of having an Empty Queue
p0 = sum(pis{1});

% Average Number of Jobs in Service given that the Queue is Empty
Ns_empty = pis{1}*Lz;
% Average Number of Jobs in Service given that the Queue is NOT Empty
Ns_no_empty = (Ns - Ns_empty);

% Average Number of Busy Servers given that the Queue is Empty
n_busy_busy_queue = (Ns_no_empty*[1 T].');
% Average Number of Busy Servers given that the Queue is NOT Empty
n_busy_no_queue = (Ns_empty*[1 T].');

% Avgerage Number of Idle Servers (Indipendently from Queue Status)
Wasted = N - B;

% Average Number of Wasted Servers (Queue is NOT Empty)
Wasted_busy_queue = N*(1-p0) - n_busy_busy_queue;
% Average Number of Wasted Servers (Queue is Empty)
Wasted_no_queue = N*p0 - n_busy_no_queue;

exeTime = toc/60;

end