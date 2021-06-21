close all; clc; clear;

n = 4; m = 1;
x = zeros(n, 1); u = 0;

P = eye(n);
p = zeros(n, 1);
costQ = eye(n);
costR = eye(m);

% Vx  = g_next_Vx( P, p, costQ, costR, x, u)
% Vxx = g_next_Vxx(P, p, costQ, costR, x, u)

for i = 1:20
p = g_next_Vx( P, p, costQ, costR, x, u)'
P = g_next_Vxx(P, p, costQ, costR, x, u)
end