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

d = sym([0; 0.00015*SymbolicEngine.v(2)]);

% T = sym([0;0.061]);
T = sym(eye(2));
SymbolicEngine.u = sym('u', [2, 1]); assume(SymbolicEngine.u, 'real');

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


% double(subs(H \ c, SymbolicEngine.q(1), 0))
% dx_0 = subs(dx,   x, zeros(size(x)));  
% dx_0 = subs(dx_0, u, zeros(size(u)));
% dx_0 = double(dx_0)

dt = 0.001;
next_x = x + dx*dt;

A = simplify(jacobian(next_x, x))
B = simplify(jacobian(next_x, u))


matlabFunction(next_x,  'File','g_f',   'Vars', {x, u});
matlabFunction(A,  'File','g_A',   'Vars', {x, u});
matlabFunction(B,  'File','g_B',   'Vars', {x, u});







