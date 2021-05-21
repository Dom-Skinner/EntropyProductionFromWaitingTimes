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

%{
P = zeros(N+1,1);
Q = zeros(N+1,1);
R = zeros(N+1,1);

for l = -N/2:N/2
    P(l+offset) = sum( exp(2*pi*1j*l*x).* b .* log(b./a))/length(x);
    Q(l+offset) = sum( exp(2*pi*1j*l*x).* (-2*p.*b./a + (log(a./b) + 1) ) )/length(x);
    R(l+offset) = sum( exp(2*pi*1j*l*x).* (2*p.*(log(b./a) + 1) - a./b))/length(x);
end
%}
%%{
make_full = @(fk) [conj(fk(end:-1:1)); fk(2:end)];

P = ifft(b .* log(b./a));
P = make_full(conj(P(1:offset)'));

Q = ifft(-2*p.*b./a + (log(a./b) + 1) );
Q = make_full(conj(Q(1:offset)'));

R = ifft(2*p.*(log(b./a) + 1) - a./b);
R = make_full(conj(R(1:offset)'));
%}
end