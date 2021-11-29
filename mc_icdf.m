function sample_mean = mc_icdf(n,alpha,beta);

syms t x;
%def PDF
c = gamma(alpha + beta)/(gamma(alpha)*gamma(beta));
f = c*t^(alpha-1)*(1-t)^(beta-1); %pdf wrt t

Us = rand(n,1); %size n array of random numbers from U(0,1)
vals = zeros(n,1); %the n uniform samples on the cdf


%% Newton to solve each x

for i=1:n
    %want to solve int_0^x f(t)dt = u for x
    eqn = int(f,t,0,x) == Us(i);
    vals(i) = solve(eqn,x);
end


%vals are the n uniform samples on the cdf

V = int(1,t,0,x);
sample_mean = (V/n)*sum(vals);

end
