function T = MatrixJacobian(M, x)

[n1, n2] = size(M);
n3 = length(x);

T = sym(zeros(n1, n2, n3));
t = jacobian(M(:), x);
for i = 1:n3
    T(:, :, i) = simplify(reshape(t(:, i), [n1, n2]));
end

end