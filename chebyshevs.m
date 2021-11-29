%%Approximate with first 6 chebyshev polynomials on [-1,1] at t

function f_final = chebyshevs()
n = 5; %highest degree chebyshev is 5

syms x;
f = exp(x);
%%k=0 to 5 chebyshev polynomials, initializations
q = chebyshevT([0, 1, 2, 3, 4, 5], x);
c = zeros(1,6);
l = n+1;

%get roots of polynomial
x_j = zeros(1,6); %preallocate
for j=0:n
    x_j(j+1) = cos((j+.5)*pi/(n+1));
end

%Calc c_k using the roots
a = 2/(n+1);
for i=1:l
    Ts = subs(q(i), x, x_j); %T_k(x_j)
    fs = subs(f, x, x_j);
    cs = times(fs,Ts); %elementwise mult to get f(x_j)T_k(x_j) %%FIXTHIS
    c(i) = a*sum(cs); %sums of f(x_j)T_k(x_j) times constant
end

%for only degree 0, the coefficient is supposed to be halved from the formula
c(1) = c(1)/2;

%% Alg: step 1
f_a2 = c(6); %f_5 = c_5

%define r and t
r = 2;
t = 1;
%since s is 0 for the chebyshev polynomials on [-1,1], it is not included

%step 2: get f_4 and f_5
f_a1 = c(5) + f_a2*(r*x); 

for i=2:4 %get f(i) from f(i+1) + f(i+2)
    f_a = c(6-i) + f_a1*(r*x) - f_a2*t;

    %update recursive expressions
    f_a2 = f_a1;
    f_a1 = f_a; 
end

f_a = c(1) + f_a1*(x) - f_a2*t;

f_final = f_a;

%% (b): graph them 
% value = subs(f_final, x, t);

xpoints = linspace(-1,1,200);
plot(xpoints, exp(xpoints));
hold on;
plot(xpoints, subs(f_final, x, xpoints));
legend("e^t", "Chebyshev approx");
title("Chebyshev Approximation of e^t");

%% (c): determine error at t=0

f_approx_0 = subs(f_final, x, 0);
f_0 = subs(f, x, 0);
%calc abs error
error = double(abs(f_approx_0 - f_0));
X = ['Absolute error is ', num2str(error), '.'];
disp(X)

%% (d): determine integrted squared error

f_err = (f_approx_0 - f)^2;
int_err = double(int(f_err, x, -1, 1));
Y = ['Integrated squared error is ', num2str(int_err), '.'];
disp(Y)

end


