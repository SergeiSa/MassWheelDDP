function M = VectorTensorMultiply(v, T)

n1 = size(T, 1);
n2 = size(T, 2);
n3 = size(T, 3);
M = sym(zeros(n2, n3));

if length(v) ~= n1
    error('dimentions inconsistent')
end

for i = 1:n3
    M(:, i) = simplify( v*T(:, :, i) );
end

end