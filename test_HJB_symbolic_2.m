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
p = sym('p', [1, n]); assume(p, 'real');

costQ = sym('Q', [n, n]); assume(costQ, 'real');
costR = sym('R', [m, m]); assume(costR, 'real');

HJB_Q = x' * costQ * x + u' * costR * u + ...
    next_x' * P * next_x + p * next_x;

HJB_Q = simplify(HJB_Q);
% vpa(HJB_Q, 3)

Qx  = simplify(jacobian(HJB_Q, x));
Qu  = simplify(jacobian(HJB_Q, u));
Quu = simplify(jacobian(Qu,    u));
Qux = simplify(jacobian(Qu,    x));
Qxx = simplify(jacobian(Qx,    x));


matlabFunction(Qx,  'File','g_Qx',   'Vars', {P, p, costQ, costR, x, u});
matlabFunction(Qu,  'File','g_Qu',   'Vars', {P, p, costQ, costR, x, u});
matlabFunction(Quu, 'File','g_Quu',  'Vars', {P, p, costQ, costR, x, u});
matlabFunction(Qux, 'File','g_Qux',  'Vars', {P, p, costQ, costR, x, u});
matlabFunction(Qxx, 'File','g_Qxx',  'Vars', {P, p, costQ, costR, x, u});
























