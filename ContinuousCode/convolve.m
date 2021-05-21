function R = convolve(fk,gk)
%this function represents the integral of f(x)*g(x)dx over the unit
%interval
ffull = [conj(fk(end:-1:1)); fk(2:end)];
gfull = [conj(gk(end:-1:1)); gk(2:end)];
R = real(sum ( ffull .* gfull(end:-1:1)));
end