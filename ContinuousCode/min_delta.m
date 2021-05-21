function [x_f,fval,exitflag] = min_delta(t2,x1_val,ctol,otol) 
syms l x0 x1 b0 b1 F C0 C1 C2 C3 a0 a1

assert(isscalar(t2),'Error: Non-scalar input')
assert(t2 > 0)

x0_val = 0; % WLOG
% x = [l,F,b0,b1] all positive
nvars = 4;
UB = (1e8)*ones(nvars,1);
LB = zeros(nvars,1);
LB(end-1:end) = 0;


%% Solve symbolically for a and p
conds = [C0 + C1 *exp( l*x0) == C2 + C3*exp(l*x0), ...
        C0 + C1*exp(l*(-1 + x1)) == C2 + C3*exp(l*x1), ...
        (a0* l)/2 - C1* exp(l *x0) *l + C3 *exp(l *x0)*l - b0*(C2 + C3* exp(l* x0))* l == 0, ...
        (a1*l)/2 + C1* exp(l *(-1 + x1)) *l - C3 *exp(l* x1)* l -  b1 * (C2 + C3* exp(l *x1))*l == 0, ...
        a0 + a1  ==1/F, ...
        (C1 *(exp(l *x0) - exp(l *(-1 + x1))))/l + (C3 *(-exp(l*x0) + exp(l*x1)))/l + C0*(1 + x0 - x1) + C2*(-x0 + x1) == 1/2];
p_fun = solve(conds,[C0,C1,C2,C3,a0,a1]);
a0_val = @(x) double(subs(p_fun.a0,[l,x0,x1,b0,b1,F],x));
a1_val = @(x) double(subs(p_fun.a1,[l,x0,x1,b0,b1,F],x));
p0_val  = @(x) double(subs(p_fun.C2 + p_fun.C3*exp(l*x0),[l,x0,x1,b0,b1,F],x));
p1_val  = @(x) double(subs(p_fun.C2 + p_fun.C3*exp(l*x1),[l,x0,x1,b0,b1,F],x));
%% solve symbolically for t2
syms C4 C5 C6 C7
conds = [C4+C5*exp(-l* x0)-x0/l==C6+C7*exp(-l* x0)-x0/l,...
    C4+C5*exp(-l*(-1+x1))-(-1+x1)/l==C6+C7*exp(-l*x1)-x1/l,...
    C5*exp(-l*x0)*l-C7*exp(-l*x0)*l-b0*l*(C6+C7*exp(-l*x0)-x0/l)==0,...
    0==-C5*exp(-l*(-1+x1))*l+C7*exp(-l*x1)*l-b1*l*(C6+C7*exp(-l *x1)-x1/l)];
Y = solve(conds,[C4,C5,C6,C7]);

C4 = Y.C4;
C5 = Y.C5;
C6 = Y.C6;
C7 = Y.C7;

syms C8 C9 C10 C11
conds2 = [C8+C9*exp(-l*x0)+(2*C5*exp(-l*x0)*x0)/l-((2*C4+2/l^2)*x0)/l+x0^2/l^2==...
    C10+C11*exp(-l*x0)+(2*C7*exp(-l*x0)*x0)/l-((2*C6+2/l^2)*x0)/l+x0^2/l^2,...
C8+C9*exp(-l*(-1+x1))+(2*C5*exp(-l*(-1+x1))*(-1+x1))/l-((2*C4+2/l^2)*(-1+x1))/l+(-1+x1)^2/l^2==...
 C10+C11*exp(-l*x1)+(2*C7*exp(-l*x1)*x1)/l-((2*C6+2/l^2)*x1)/l+x1^2/l^2,...
-((2*C5*exp(-l*x0))/l)+(2*C7*exp(-l*x0))/l+(2*C4+2/l^2)/l-(2*C6+2/l^2)/l-...
C11*exp(-l*x0)*l+C9*exp(-l*x0)*l+2*C5*exp(-l*x0)*x0-2*C7*exp(-l*x0)*x0-...
b0*l*(C10+C11*exp(-l*x0)+(2*C7*exp(-l*x0)*x0)/l-((2*C6+2/l^2)*x0)/l + x0^2/l^2)==0,...
(2*C5*exp(-l*(-1+x1)))/l-(2*C7*exp(-l*x1))/l-(2*C4+2/l^2)/l+(2*C6+2/l^2)/l-C9*exp(-l*(-1+x1))*l+...
C11*exp(-l*x1)*l-2*C5*exp(-l*(-1+x1))*(-1+x1)+(2*(-1+x1))/l^2+2*C7*exp(-l*x1)*x1-...
(2*x1)/l^2-b1*l*(C10+C11*exp(-l*x1)+(2*C7*exp(-l*x1)*x1)/l-((2*C6+2/l^2)*x1)/l+x1^2/l^2)==0];

Y = solve(conds2,[C8,C9,C10,C11]);

t2_val0 = @(x) double(subs(Y.C10+Y.C11*exp(-l*x0)+(2*C7*exp(-l*x0)*x0)/l-((2*C6+2/l^2)*x0)/l+x0^2/l^2,[l,x0,x1,b0,b1,F],x));
t2_val1 = @(x) double(subs(Y.C10+Y.C11*exp(-l*x1)+(2*C7*exp(-l*x1)*x1)/l-((2*C6+2/l^2)*x1)/l+x1^2/l^2,[l,x0,x1,b0,b1,F],x));

t2_diff = @(l_val,F_val,b0_val,b1_val,x0_val,x1_val) l_val^2/F_val*(a0_val([l_val,x0_val,x1_val,b0_val,b1_val,F_val])*t2_val0([l_val,x0_val,x1_val,b0_val,b1_val,F_val])+...
    a1_val([l_val,x0_val,x1_val,b0_val,b1_val,F_val])*t2_val1([l_val,x0_val,x1_val,b0_val,b1_val,F_val]))-t2;
%% set up necessary initial conditions
sig = @(l_val,F_val,b0_val,b1_val,x0_val,x1_val) F_val*(l_val/2+p0_val([l_val,x0_val,x1_val,b0_val,b1_val,F_val])*b0_val*log(b0_val/a0_val([l_val,x0_val,x1_val,b0_val,b1_val,F_val]))+p1_val([l_val,x0_val,x1_val,b0_val,b1_val,F_val])*b1_val*log(b1_val/a1_val([l_val,x0_val,x1_val,b0_val,b1_val,F_val]))+...
        (a0_val([l_val,x0_val,x1_val,b0_val,b1_val,F_val])*log(a0_val([l_val,x0_val,x1_val,b0_val,b1_val,F_val])/b0_val)+a1_val([l_val,x0_val,x1_val,b0_val,b1_val,F_val])*log(a1_val([l_val,x0_val,x1_val,b0_val,b1_val,F_val])/b1_val))/2);
ObjectiveFunction = @(x) sig(x(1),x(2),x(3),x(4),x0_val,x1_val);

ConstraintFunction = @(x) deal([],t2_diff(x(1),x(2),x(3),x(4),x0_val,x1_val)); % the one non-linear constraint

% Set up the problem with an initial guess and solve
l_init = 6;
F_inir = 1.25;
b0_init = 0.15;
b1_init = 0.8;
x_init = [l_init;F_inir;b0_init;b1_init];

options = optimoptions('fmincon','Display','iter','ConstraintTolerance',ctol,...
    'OptimalityTolerance',otol,'MaxFunctionEvaluations', 1e+04);
%
%[x,fval,exitflag,~]  = fmincon(ObjectiveFunction,x0,A,b,[],[],LB,UB,ConstraintFunction,options);

problem = createOptimProblem('fmincon','x0',x_init,'objective',ObjectiveFunction,...
        'lb',LB,'ub',UB, 'nonlcon',ConstraintFunction,'options',options);
%gs = GlobalSearch('MaxTime',3600);
ms = MultiStart;
[x_f,fval,exitflag,~] = run(ms,problem,4);

end