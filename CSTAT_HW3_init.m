%% CStats HW init file
syms x

%% 3.1 
chebyshevs(); %running the file chebyshevs() gives all desired outputs

%% 3.2 approx function two diff ways, plot them. 

ns = [1, 2, 4, 8, 16]; %set of ns to do this for

%create a subgraph using these with the errors printed
for i=1:len(ns)
    n = ns(i);
    [plots(i), e1(i), e2(i)] = approx(n);
end

%create subplot and print errors

%% 3.3 use trapezoid rule to approx the function, print table for given x0s.

x0s = linspace(0,8,9);
f = @(x)(2^x);
x0s = arrayfun(f,x0s);

fx = zeros(9,1);
ns = zeros(9,1);

for i=1:8
    [fx(i), ns(i)] = trapezoid_acc(x0s(i));
end

A = [x0s.', fx, ns];
varTypes = {'int8', 'double', 'int8'};
varNames = ["x", "J_0(x)", "Minimum n"];
a2t = array2table(A,"VariableNames",varNames);
disp(a2t);

%% 3.4

%(a) compute these values with 8th order adaptive legendre-gauss quad.
ans = adaptive_quad(.25,3,3);
A = ['F(.25,3,3) = ', num2str(ans), '.'];
disp(A)
ans = adaptive_quad(.5,4,5);
A = ['F(.5,4,5) = ', num2str(ans), '.'];
disp(A)

% do the same thing with montecarlo cdf method
ns = [10, 100, 1000, 10000, 100000];

for n=len(ns)
    n = ns(i);
    sample_mean = mc_icdf(n,5,6);
    Value of sample mean at N = 
    A = ['Value of sample mean at N = ', num2str(n), ' is ', sample_mean, '.'];
    disp(A)
end
real = 5/(5+6);
A = ['We see that convergence to real mean of ', num2str(real), ' is approximately 1/sqrt(N).'];
    disp(A)

