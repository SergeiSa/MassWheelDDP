clear; 
LinkArray = SRD_get('LinkArray');

g0 = sym('g');    assume(g0, 'real');
d0 = sym('d');    assume(d0, 'real');
tau = sym('tau'); assume(tau, 'real');
L = sym('L'); assume(L, 'real');
u = sym('u'); assume(u, 'real');
m = sym('m', [2, 1]); assume(m, 'real');
I = sym('I', [2, 1]); assume(I, 'real');

LinkArray(2).Mass = m(1);
LinkArray(3).Mass = m(2);
LinkArray(2).Inertia = eye(3)*I(1);
LinkArray(3).Inertia = eye(3)*I(2);
LinkArray(2).RelativeCoM = [0;0;L/2];
LinkArray(2).RelativeFollower = [0;0;L];

SymbolicEngine = SRDSymbolicEngine('LinkArray', LinkArray, 'Casadi', false);
SymbolicEngine.InitializeLinkArray();

SRD_dynamics_derive_JacobiansForLinkArray('SymbolicEngine', SymbolicEngine);

H = SRD_dynamics_derive_JSIM('SymbolicEngine', SymbolicEngine);

[in, dH] = SRD_dynamics_derive_GeneralizedInertialForces_via_dH(...
    'SymbolicEngine', SymbolicEngine, ...
    'JointSpaceInertiaMatrix', H);
g = SRD_dynamics_derive_GeneralizedGravitationalForces(...
    'SymbolicEngine', SymbolicEngine, ...
    'GravitationalConstant', [0; 0; -g0]);


d = sym([0;0;0;d0]);

T = sym([0;0;0;tau]);
c = in + g + d;

f = simplify( H \ (T*u - c) )
Aq = simplify(jacobian(f, SymbolicEngine.q))
Av = simplify(jacobian(f, SymbolicEngine.v))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
constraint = [SymbolicEngine.q(1); SymbolicEngine.q(2)];

k = length(constraint);
n = length(SymbolicEngine.q);

F = jacobian(constraint, SymbolicEngine.q);
dF = zeros(size(F));

M = [H, -F';
     F,  zeros(k,k)];
RHS = [T*u - c; -dF * SymbolicEngine.v];
map = [eye(n), zeros(n, k)];

fc = simplify( map*(M \ RHS) )
Aqc = simplify(jacobian(fc, SymbolicEngine.q))
Avc = simplify(jacobian(fc, SymbolicEngine.v))

