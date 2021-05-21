function sig = entropy_rate(a,c,p,F,D)
x = linspace(0,1,400);
x = x(1:end-1);

rs = @(fk,w) 2*real(sum( fk(2:end) .* w)) + real(fk(1));
w = exp( 1j * 2* pi* (1:length(a)-1)' * x);
%rs_p = @(fk) 2*real(sum( (1j * 2* pi* (1:length(fk)-1)').*fk(2:end) .* exp( 1j * 2* pi* (1:length(fk)-1)' * x)));
%rs_p2 = @(fk) 2*real(sum( (1j * 2* pi* (1:length(fk)-1)').^2  .*fk(2:end) .* exp( 1j * 2* pi* (1:length(fk)-1)' * x)));
a = rs(a,w);
c = rs(c,w);
p = rs(p,w);
%Fp = rs_p(F);
F = rs(F,w);
%Dp = rs_p(D);
%Dpp = rs_p2(D);
D = rs(D,w);



%sig = sum( 0.5*a.*log(a./c) + p.*c.*log(c./a) + p.*( F.^2 ./D + Fp - 3*F.*Dp./D + Dp.^2 ./D - Dpp))/length(x);
sig = sum( 0.5*a.*log(a./c) + p.*c.*log(c./a) + 0.5.*( F.^2 ./D) )/length(x);
end