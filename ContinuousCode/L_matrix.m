function L = L_matrix(Fk,Dk,bk)
N = 2*length(Fk)-2;
offset = length(Fk);
%{
% For general matrix with F, D non constant
L = zeros(N+1,N+1);
for r = -(N/2):(N/2)
    for k = -(N/2):(N/2)
            if (r-k + offset > 0) && (r-k +offset <= N+1)
                if r-k < 0
                    L(r + offset,k+offset) = conj(Fk(abs(r-k)+1))*(2*pi*1j*k) ...
                        - conj(Dk(abs(r-k)+1))*(2*pi*k)^2-conj(ck(abs(r-k) +1));
                else
                    L(r + offset,k+offset) = Fk(r-k+1)*(2*pi*1j*k) - Dk(r-k+1)*(2*pi*k)^2-ck(r-k +1);
                end
            end
    end
end
%}

% This faster version only works for F, D const.
make_full = @(fk) [conj(fk(end:-1:1)); fk(2:end)];
b_conj= make_full(bk)';
L = diag(2*pi*Fk(1)*1j*(-N/2:N/2) -Dk(1)*(2*pi)^2 *(-N/2:N/2).^2);
offset = N/2+1;
for k = -(N/2):(N/2)
    l = -N/2:N/2;
    dxs = l((l - k >= -N/2) & (l-k <= N/2)) + offset;   
    L(k+offset,dxs) = L(k+offset,dxs) - b_conj( dxs-k);
end

end