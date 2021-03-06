close all; clc; clear;


n = 4; m = 2;
x = zeros(n, 1); u = 0;

P  = eye(n);
costQ = eye(n);
costR = eye(m);

A = g_A(x, u);
B = g_B(x, u);

Count = 10000;
change_log = zeros(Count, 1);
K_log = zeros(m, n, Count);

for i = 1:Count
    K = pinv(costR + B'*P*B) * B' * P * A;
    P_next = costQ + K' * costR * K + (A-B*K)' * P* (A-B*K);
%     K = pinv(costR + B'*P*B) *B'*P*A;
%     P_next = costQ + A'*P*(A - B*K);
    
    change_log(i) = norm(P - P_next);
    K_log(:, :, i) = K;
    
    P = P_next;
end

figure()
subplot(1, 2, 1); semilogy(change_log)
K
P
eig(A - B*K)

x_log = zeros(Count, n);
x = 0.01*randn(n, 1);
for i = 1:Count
    x_log(i, :) = x;
    K = K_log(:, :, Count - i + 1);
    x = g_f(x, -K*x);
end
subplot(1, 2, 2)
plot(x_log)

%dt = 0.001;
% odefnc = @(t, x) g_f(x, K*x);
% x0 = 0.01*randn(n, 1);
% tf = Count * dt;
% [TOUT, XOUT] = ode45(odefnc,[0, tf],x0);
% plot(TOUT, XOUT)
