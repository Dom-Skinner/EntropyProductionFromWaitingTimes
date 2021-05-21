function sig = entropy_rate(a,c,p,F,D)
x = linspace(0,1,400);
x = x(1:end-1);

rs = @(fk,w) 2*real(sum( fk(2:end) .* w)) + real(fk(1));
w = exp( 1j * 2* pi* (1:length(a)-1)' * x);
a = rs(a,w);
c = rs(c,w);
p = rs(p,w);
F = rs(F,w);

D = rs(D,w);

sig = sum( 0.5*a.*log(a./c) + p.*c.*log(c./a) + 0.5.*( F.^2 ./D) )/length(x);
end