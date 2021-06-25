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

dt = 0.001;
next_x = x + dx*dt;


%%%%%%%%%%%%%%%%%%%%%%%
   
  
Vxx = sym('P', [n, n]); assume(Vxx, 'real');
Vx  = sym('p', [1, n]); assume(Vx, 'real');

costQ = sym('Q', [n, n]); assume(costQ, 'real');
costR = sym('R', [m, m]); assume(costR, 'real');

lx = x'*costQ;
lu = u'*costR;

lxx = costQ;
luu = costR;
lxu = zeros([size(costQ, 1), size(costR, 1)]);

fx = simplify(jacobian(next_x, x));
fu = simplify(jacobian(next_x, u));
fxx = MatrixJacobian(fx, x);
fxu = MatrixJacobian(fx, u);
fux = MatrixJacobian(fu, x);
fuu = MatrixJacobian(fu, u);

Qx = simplify(lx + Vx*fx) %multiplication dimentions?
Qu = simplify(lu + Vx*fu) %multiplication dimentions?
Qxx = simplify(lxx + fx' * Vxx * fx + VectorTensorMultiply(Vx, fxx))
Qxu = simplify(lxu + fx' * Vxx * fu + VectorTensorMultiply(Vx, fxu))
Quu = simplify(luu + fu' * Vxx * fu + VectorTensorMultiply(Vx, fuu))


matlabFunction(Qx,  'File','g2_Qx',   'Vars', {Vxx, Vx, costQ, costR, x, u});
matlabFunction(Qu,  'File','g2_Qu',   'Vars', {Vxx, Vx, costQ, costR, x, u});
matlabFunction(Quu, 'File','g2_Quu',  'Vars', {Vxx, Vx, costQ, costR, x, u});
matlabFunction(Qxu, 'File','g2_Qxu',  'Vars', {Vxx, Vx, costQ, costR, x, u});
matlabFunction(Qxx, 'File','g2_Qxx',  'Vars', {Vxx, Vx, costQ, costR, x, u});







