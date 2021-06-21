close all; clc; clear;

n = 4; m = 1;
x = zeros(n, 1); u = 0;

Vxx  = eye(n);
Vx   = zeros(1, n);
costQ = eye(n);
costR = eye(m);


for i = 1:20
Qx = g_Qx(Vxx, Vx, costQ, costR, x, u);
Qu = g_Qu(Vxx, Vx, costQ, costR, x, u);
Quu = g_Quu(Vxx, Vx, costQ, costR, x, u);
Qux = g_Qux(Vxx, Vx, costQ, costR, x, u);
Qxx = g_Qxx(Vxx, Vx, costQ, costR, x, u);

Vx = Qx - (Qu' / Quu) * Qux;
Vxx = Qxx - (Qux' / Quu) * Qux;
eig(Vxx)
end

