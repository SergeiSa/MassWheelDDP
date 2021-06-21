function [N_table, G_table, An_table, Bn_table, cn_table, xn_table] ... 
    = SRD_ConstrainedLinearModel_EstimatorController_v1(varargin)
Parser = inputParser;
Parser.FunctionName = 'SRD_LinearModel_GenerateTable';
Parser.addOptional('Handler_dynamics_generalized_coordinates_model', []);
Parser.addOptional('Handler_dynamics_Linearized_Model', []);
Parser.addOptional('Handler_Constraints_Model', []);
Parser.addOptional('x_table', []);
Parser.addOptional('C', []);
Parser.addOptional('tol', 10^(-5));

Parser.parse(varargin{:});


Count = size(Parser.Results.x_table, 2);
n = Parser.Results.Handler_dynamics_generalized_coordinates_model.dof_configuration_space_robot;
m = Parser.Results.Handler_dynamics_generalized_coordinates_model.dof_control;
k = Parser.Results.Handler_Constraints_Model.dof_Constraint;
    
C   = Parser.Results.C;
tol = Parser.Results.tol;

N_table = zeros(n, nn, Count);
G_table = zeros(2*k, n, Count);
An_table = zeros(nn, nn, Count);
Bn_table = zeros(nn, m, Count);
cn_table = zeros(nn, Count);
xn_table = zeros(nn, Count);

for i = 1:Count
    
    x = Parser.Results.x_table(:, i);
    q = x(1:n);
    v = x((n + 1):end);
    u = zeros(m, 1);
    
    H = Handler_dynamics_generalized_coordinates_model.get_joint_space_inertia_matrix(q);
    c = Handler_dynamics_generalized_coordinates_model.get_bias_vector(q, v);
    T = Handler_dynamics_generalized_coordinates_model.get_control_map();
    iH = pinv(H);
    f0 = [v; iH*(T*u - c)];
    
    A = Handler_dynamics_Linearized_Model.get_A(q, v, u, iH);
    B = Handler_dynamics_Linearized_Model.get_B(q, v,    iH);
    g = f0 - A * x - B * u;
    
    F  = Handler_Constraints_Model.get_Jacobian(q);
    dF = Handler_Constraints_Model.get_Jacobian_derivative(q);
    G = [zeros(k, n), F; F, dF];
    G = svd_suit(G, tol);    
    N = G.null;
    R = G.row_space;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %D1*dx/dt + D2*x = 0
    D1 = [eye(n), zeros(n)];
    D2 = [zeros(n), eye(n)];

    GG = [G.self, zeros(2*k, 2*n); D1, D2];
    GG = svd_suit(GG, tol);

    V = GG.null; V = V((2*n+1):end, :);

    R_ES = svd_suit( ( pinv(N'*A)*(N'*A) * V*pinv(V) * R), tol);
    R_ES = R_ES.orth;
    
    E = [N, R_ES];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nn_1 = size(N,  2);
    nn_2 = size(R_ES, 2);
    nn = nn_1+nn_2;
    mm = size(C, 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    desired.M = svd_suit([N'*A*N, N'*B]);
    desired.N = desired.M.null;
    desired.map.z = [eye(nn_1),   zeros(nn_1, m)];
    desired.map.u = [zeros(m, nn_1), eye(m)];
    desired.z_particular = desired.map.z * desired.M.pinv * N'*g;
    desired.Projector = (desired.map.z*desired.N) * pinv(desired.map.z*desired.N);
    desired.z_desired = desired.z_particular + desired.Projector * (desired.z_desired_0 - desired.z_particular);

    
    G_table(:, :, i) = G;
    N_table(:, :, i) = N;
    
    An_table(:, :, i) = N' * Parser.Results.A_table(:, :, i) * N;
    Bn_table(:, :, i) = N' * Parser.Results.B_table(:, :, i);
    cn_table(:, i)    = N' * Parser.Results.c_table(:, i);
    
    xn_table(:, i) = N' * Parser.Results.x_table(:, i);
    dxn_table(:, i) = N' * Parser.Results.dx_table(:, i);
    
end

end