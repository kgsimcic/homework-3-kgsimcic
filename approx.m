function [e1, e2] = approx(n) %approximate the function according to n points

syms t x;

f = cos(x);

%% (a) polynomial interpolation
xs = linspace(0,pi,n);
y = cos(xs);

p = polyfit(xs,y,n-1); %p gives p_1x^n-1 + .... + p_n+1

%assemble actual polynomial to evaluate as p_x
p_x = 0;
for i=0:n-1
    p_x = p(i+1)*x^(n-1-i) + p_x; %creates actual polynomial
end

%% (b) project onto set of polies P_0, ..., P_n-1 -- scaled/trans legendre polies

%x is a function of t
f = @(t)(cos((pi/2)*(t+1)));

for i=0:n-1
    q(i+1) = legendreP(i,t); %evals at translated points
    num = (2*i+1)/2;
    c(i+1) = num*int(q(i+1)*f, t, -1, 1);
end

%assemble polynomial = sum_0^(n-1) c(i+1)q(i+1)

L_x = 0;
for i=0:n-1
    L_x = c(i+1)*q(i+1) + L_x; %creates actual polynomial sum
end

%% errors

L_x = subs(L_x, t, (2/pi)*x - 1); %translates L_x back

%change f back into terms of x
f = cos(x);

Q = ['With n = ', num2str(n), ','];
disp(Q)

fxn1 = (p_x - f)^2;
e1 = double(sqrt(int(fxn1, x, 0, pi)));
A = ['L2 error for (a) was ', num2str(e1)];
disp(A)

fxn2 = (L_x - f)^2;
e2 = double(sqrt(int(fxn2, x, 0, pi)));
B = ['L2 error for (b) was ', num2str(e2)];
disp(B)

end