function [P,Q,R] = entropy_grad(a,b,p)
N = 2*length(a)-2;
offset = N/2+1;


x = linspace(0,1,400);
x = x(1:end-1);

rs = @(fk,w) 2*real(sum( fk(2:end) .* w)) + real(fk(1));
w = exp( 1j * 2* pi* (1:length(a)-1)' * x);
a = rs(a,w);
b = rs(b,w);
p = rs(p,w);


make_full = @(fk) [conj(fk(end:-1:1)); fk(2:end)];

P = ifft(b .* log(b./a));
P = make_full(conj(P(1:offset)'));

Q = ifft(-2*p.*b./a + (log(a./b) + 1) );
Q = make_full(conj(Q(1:offset)'));

R = ifft(2*p.*(log(b./a) + 1) - a./b);
R = make_full(conj(R(1:offset)'));
end