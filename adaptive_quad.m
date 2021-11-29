%%Want approx CDF of below dist with inputs alpha, beta, t:
%NOTE: I switched t and x because I use syms x in matlab. t is input.

function val = adaptive_quad(x,alpha,beta)

syms t;
%def PDF
c = gamma(alpha + beta)/(gamma(alpha)*gamma(beta));
f = c*t^(alpha-1)*(1-t)^(beta-1); %pdf wrt t

% get real integral
f_real = int(f,t,0,x);
a = 0;
b = x;

%perform initial Gauss_legendre quad: Used lgwt from the internet. 
[nodes, weights] = lgwt(8,a,b); %8 cause 8th order
fs = subs(f, t, nodes); %evaluate function at the nodes first
Q = double(sum(fs.*weights)); %dotting this with the weights and summing gives the approximation

tau = abs(Q - f_real);

%want to approx int_0^t f(x) dx with precision of 10^(-8)
epsilon = 10^(-8);

while epsilon < tau
    m = (a+b)/2;

    [nodes1, weights1] = lgwt(8,a,m); [nodes2, weights2] = lgwt(8,m,b);
    fs1 = f(nodes1); fs2 = f(nodes2);
    Q = sum(fs1.*weights1) + sum(fs2.*weights2);
    Q = double(Q);
end

val = Q;
end









