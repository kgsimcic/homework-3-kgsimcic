%calculates required n for min accuracy of 10^-12 of approximating
%1/pi * int_0^pi cos(x_0cos(theta))dtheta by using composite trap rule

function [final_f, n] = trapezoid_acc(x_0)

n = 1;
x = sym('x');
f = @(x)(cos(x_0*cos(x)));
%f = cos(x_0*cos(x));

f_real = (1/pi)*int(f, x, 0, pi); %calc real value to check
Fvpa = double(vpa(f_real));
%make first approximation with n=1:
f_approx = (1/2)*(f(pi) + f(0)); %cancelled 1/pi and pi/2

while abs(vpa(f_approx - f_real)) > 10^(-12)
    %increment n and h
    n = n + 1;
    h = pi/n;

    %sum f(x_i):
    fs = zeros(n+1,1);
    xs = linspace(0,pi, n+1); %gives middle points
    fs = (1/pi)*f(xs);
    fs = sum(fs);

    %new approx for n+1
    f_approx = h*(fs) - (h/pi)*((f(pi) + f(0))/2);
end

final_f = f_approx;

end