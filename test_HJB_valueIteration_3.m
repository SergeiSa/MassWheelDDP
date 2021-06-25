close all; clc; clear;

n = 4; m = 2;
x = zeros(n, 1); u = zeros(m, 1);

Vxx  = eye(n);
Vx   = zeros(1, n);
costQ = eye(n);
costR = eye(m);


Count = 2000;
change_log = zeros(Count, 2);
K_log = zeros(m, n, Count);
k_log = zeros(m,    Count);

for i = 1:Count
Qx  = g2_Qx (Vxx, Vx, costQ, costR, x, u);
Qu  = g2_Qu (Vxx, Vx, costQ, costR, x, u);
Quu = g2_Quu(Vxx, Vx, costQ, costR, x, u);
Qux = g2_Qxu(Vxx, Vx, costQ, costR, x, u);
Qxx = g2_Qxx(Vxx, Vx, costQ, costR, x, u);

Vx_next = Qx - Qu * pinv(Quu) * Qux';
Vxx_next = Qxx - Qux * pinv(Quu) * Qux';

change_log(i, :) = [norm(Vxx_next - Vxx), norm(Vx_next - Vx)];

K_log(:, :, i) = pinv(Quu) * Qux';
k_log(:,    i) = pinv(Quu) * Qu';

Vx = Vx_next;
Vxx = Vxx_next;
end

figure()
subplot(1, 2, 1); semilogy(change_log)


x_log = zeros(Count, n);
x = 0.01*randn(n, 1);
for i = 1:Count
    x_log(i, :) = x;
    
    K = K_log(:, :, Count - i + 1);
    k = k_log(:,    Count - i + 1);
    
    x = g_f(x, -K*x - k);
end
subplot(1, 2, 2)
plot(x_log)