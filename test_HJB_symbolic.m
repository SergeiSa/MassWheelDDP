close all; clear; %clear classes;
clc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics

LinkArray = SRD_get('LinkArray');

SymbolicEngine = SRDSymbolicEngine('LinkArray', LinkArray, 'Casadi', false);
SymbolicEngine.InitializeLinkArray();

SRD_dynamics_derive_JacobiansForLinkArray('SymbolicEngine', SymbolicEngine);

H = SRD_dynamics_derive_JSIM('SymbolicEngine', SymbolicEngine);

[in, dH] = SRD_dynamics_derive_GeneralizedInertialForces_via_dH(...
    'SymbolicEngine', SymbolicEngine, ...
    'JointSpaceInertiaMatrix', H);

g = SRD_dynamics_derive_GeneralizedGravitationalForces(...
    'SymbolicEngine', SymbolicEngine, ...
    'GravitationalConstant', [0; 0; -9.8]);

d = sym([0;0.00015]);

T = sym([0;0.061]);

c = in + g + d;

description = SRD_generate_dynamics_generalized_coordinates_model(...
    'SymbolicEngine', SymbolicEngine, ...
    'H', H, ...
    'c', c, ...
    'T', T, ...
    'Symbolic_ToOptimizeFunctions', true, ...
    'Casadi_cfile_name', 'g_dynamics_generalized_coordinates', ...
    'FunctionName_H', 'g_dynamics_H', ...
    'FunctionName_c', 'g_dynamics_c', ...
    'FunctionName_T', 'g_dynamics_T', ...
    'Path', 'Dynamics/');

Handler_dynamics_generalized_coordinates_model = SRD_get_handler__dynamics_generalized_coordinates_model('description', description);
SRD_save(Handler_dynamics_generalized_coordinates_model, 'Handler_dynamics_generalized_coordinates_model');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dof = Handler_dynamics_generalized_coordinates_model.dof_configuration_space_robot;

x = sym('x', [dof * 2, 1]);  assume(x, 'real');
u = SymbolicEngine.u;        assume(u, 'real');

n = length(x);  
m = length(u); 

dx = [SymbolicEngine.v;
      H \ (T*u - c)];
  
dx = subs(dx, SymbolicEngine.q, x(1:dof));  
dx = subs(dx, SymbolicEngine.v, x((dof+1):end));  

dt = 0.001;
next_x = x + dx*dt;


%%%%%%%%%%%%%%%%%%%%%%%
   
  
P = sym('P', [n, n]); assume(P, 'real');
p = sym('p', [n, 1]); assume(p, 'real');

costQ = sym('Q', [n, n]); assume(costQ, 'real');
costR = sym('R', [m, m]); assume(costR, 'real');

HJB_Q = x' * costQ * x + u' * costR * u + ...
    next_x' * P * next_x + p' * next_x;

HJB_Q = simplify(HJB_Q);
% vpa(HJB_Q, 3)

Qx  = simplify(jacobian(HJB_Q, x));
Qu  = simplify(jacobian(HJB_Q, u));
Quu = simplify(jacobian(Qu,    u));
Qux = simplify(jacobian(Qu,    x));
Qxx = simplify(jacobian(Qx,    x));

% vpa(Qx, 3)
% vpa(Qu, 3)
% vpa(Quu, 3)
% vpa(Qux, 3)
% vpa(Qxx, 3)

u_exp = - Quu \ (Qu + Qux * x);

M = 0.5 * [0,   Qx,   Qu;
           Qx', Qxx,  Qux';
           Qu', Qux,  Quu];

HJB_Q_new = ...
    [1; x; u_exp]' * M * [1; x; u_exp];
HJB_Q_new = simplify(HJB_Q_new);
% vpa(HJB_Q_new, 3)

p_new = simplify(jacobian(HJB_Q_new, x));
P_new = simplify(jacobian(p_new,     x));
vpa(p_new, 3)
vpa(P_new, 3)

matlabFunction(p_new, 'File','g_next_Vx',  'Vars', {P, p, costQ, costR, x, u});
matlabFunction(P_new, 'File','g_next_Vxx', 'Vars', {P, p, costQ, costR, x, u});
























