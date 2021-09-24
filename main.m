clear all; % clc;
close all;
rng(2222);

% the number of target
K = 3;
% the number of UAV
N = 16;
% the number of measurements
M = 32;
% the distance between adjacent UAVs
d = 0.5;
posUAV = [0:N-1].'*d;
% the direction of central UAV
psi = 0;


vec = @(MAT) MAT(:);
vecH = @(MAT) MAT(:).';
vec_idxV = @(VEC, idxV) VEC(idxV);
steerVec = @(ang, pos) exp(1j*2*pi*pos*sind(vecH(ang)));


% target direction range
theta_max = 45;
theta_min = -45;

% target signals
s = sqrt(1/2) * (ones(K, 1) + 1j*ones(K,1));

isFig = false;

MC_num = 10;
ang_grid = [-50:1e-3:50].';

SNR = 20;

% SNR = 30;

% 30dB, t=600
% 20dB, t=3e3
% 10dB, t=5e4
% t=10^(-0.096*SNR+5.5722)$

theta = [-30.345, 0.789, 20.456].';


noise = sqrt(1/2)*(randn(M,1)+1j*randn(M,1));


d_per = randn(N, 1)*0.1;




ang_range = [-60:0.1:60].';
for idx_test = 1:length(ang_range)
    tmp1 = steerVec(ang_range(idx_test), posUAV);
    tmp2 = steerVec(ang_range(idx_test), posUAV+d_per);
    if idx_test == 1
        AA_T = zeros(length(tmp1), length(ang_range));
        BB_T = zeros(length(tmp2), length(ang_range));
    end
    
    AA_T(:, idx_test) = tmp1;
    BB_T(:, idx_test) = tmp2;
end
est_t = pinv(kron(AA_T.',eye(N)))*vec(BB_T);
est_T = reshape(est_t, N, N)';

% measurement matrix
B = rand(M, N);
B(B>=0.5) = 1;
B(B<0.5) = -1;
C = B*diag(steerVec(psi, posUAV+d_per));
r_no = C*steerVec(theta, posUAV+d_per)*s;

AA = pinv(est_T*C');
BB = pinv(C');

t=10^(-0.096*SNR+5.5722); 

noiseVar = norm(r_no)^2/10^(SNR/10)/length(r_no);
r = r_no+ sqrt(noiseVar)* noise;

 
%% estimation DOA 
% with pertubration
cvx_solver sdpt3
cvx_begin sdp quiet
    variable G(N+1, N+1) hermitian;
    G>=0;
    G(N+1, N+1) == 1;
    trace(G) <= 1+t;
    for idx = 1:N-1
        sum(diag(G(1:N, 1:N), idx)) == 0;
    end
    minimize(norm(r-AA*G(1:N, 1+N)))
cvx_end
x = G(1:N, N+1);
% plot
sp = abs(x'*steerVec(ang_grid, posUAV));
sp = sp/max(sp);


figure; stem(theta, ones(K,1),'BaseValue', 0, 'LineWidth', 3);
hold on; plot(ang_grid, sp, 'LineWidth', 3); 

[pks,locs] = findpeaks(sp);
[~, sort_idx] = sort(pks, 'descend');
est_ang = sort(ang_grid(locs(sort_idx(1:K))),'ascend');
 
fprintf('Proposed: %4f\n', sqrt(norm(est_ang-theta)^2/K));


 
 