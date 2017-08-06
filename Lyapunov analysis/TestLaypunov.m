

%% create a nonlinear dynamic system
clc;clear
n = 6;         %% number of variable in x, x1,x2,...,xn.
d = 2;          %% degree of fi(x)
density = 0.3;


R = sprand(n,n,density);  %% generate random sparse matrix, determine the pattern of f(x).
R = full(spones(R));
rho = 0.4;
[f,x] = GenNonDynamics(R,rho,d); 


%[ prog,SOLV, SOLR] = FindLF(f,x,'CDCS')

% local region; g(x) <= 0
gamma = 0.01;
g = sum(x.^2)-gamma;

%% find Lyapunov function
Degree = 2;
% =============================================
% First, initialize the sum of squares program
prog = sosprogram(x);

% =============================================
% The Lyapunov function V(x): 
CandidateMonomails = monomials(x,[2:Degree]);         % full candidate
%CandidateMonomails = x.^2;   % diagnal case quadratic one.

[prog,V] = sospolyvar(prog,CandidateMonomails);  % diagonal case

% =============================================
% Next, define SOSP constraints
% Constraint 1 : V(x) - 0.1*(x1^2 + x2^2 + x3^2 + ... + xn^2) - 0.1*g >= 0
kappa1 = 0.01;
PolyIneq1 = V - kappa1*sum(x.^2);
prog = sosineq(prog,PolyIneq1);

% Constraint 2: -dV/dx*f - 0.1*g >= 0;
[prog,r] = sospolyvar(prog,CandidateMonomails); 
prog = sosineq(prog,r);

PolyIneq2 = r*g;
for i = 1:length(x)
    PolyIneq2 = PolyIneq2 - diff(V,x(i))*f(i);
end;
prog = sosineq(prog,PolyIneq2);




% =============================================
% And call solver
% prog = sossolve(prog);

%% call solver
% prog.solinfo.x = [];
% options.solver = 'cdcs';
% %options.solver = 'sosadmm';
%  options.params.solver = 'primal';
% prog = sossolve(prog,options);
% SOLV = sosgetsol(prog,V)
% SOLR = sosgetsol(prog,r)

%% by CDCS-SOS
prog.solinfo.x = [];
options.solver = 'cdcs';
options.params.solver = 'sos';
prog = sossolve(prog,options);
SOLV = sosgetsol(prog,V)
SOLR = sosgetsol(prog,r)

%% by SOSADMM
% prog.solinfo.x = [];
% options.solver = 'sosadmm';
% prog = sossolve(prog,options);
% SOLV = sosgetsol(prog,V)
% SOLR = sosgetsol(prog,r)

%% by SeDuMi
prog.solinfo.x = [];
options.solver = 'sedumi';
prog = sossolve(prog,options);
SOLV = sosgetsol(prog,V)
SOLR = sosgetsol(prog,r)



%% check the solution?
% 
% PolyIneq1 = SOLV - kappa1*sum(x.^2);
% SOLR = sosgetsol(prog,r);
% % 
%  PolyIneq2 = SOLR*g;
%  for i = 1:length(x)
%      PolyIneq2 = PolyIneq2 - diff(SOLV,x(i))*f(i);
%  end;
% 
%  [Q1,Z1] = findsos(PolyIneq1);
%  [Q2,Z2] = findsos(PolyIneq2);
% 
% prog1 = sosprogram(x);
% prog1 = sosineq(prog1,PolyIneq2);
% options.solver = 'cdcs';
% % options.params.Max_iter = 1000;
% % options.params.solver = 'primal';
% % options.solver = 'cdcs';
% prog1 = sossolve(prog1,options);

