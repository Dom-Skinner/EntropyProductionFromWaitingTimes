function [c,ceq,DC,DCeq] = opt_constraint(x,N,M,t2)
    make_full = @(fk) [conj(fk(end:-1:1)); fk(2:end)];
    reduce = @(fk) fk((numel(fk)-1)/2+1 :end);
    [ak,ck,Dk,Fk] = unpack(x,N,M);
    
    L = L_matrix(Fk,Dk,ck);
    ds = zeros(N/2+1,1);
    ds(1) = 1;
    T1 = (L \ make_full(ds));
    T = 2*reduce(L\ T1);
    p = -0.5*((L') \ make_full(ak));
    T2 = convolve(T,ak);
    
    u = L*make_full(T);
    v = -2 *((L') \ p);
    
    dc1_dar = 2*real(T);
    dc1_dai = 2*imag(T);
    
    dc2_dar = -real(T1);
    dc2_dai = -imag(T1);
    
    T  = make_full(T);
    dc1_dF = 2*pi*imag(sum( (-N/2:N/2)' .* (conj(v) .* u - 2*conj(p).*T)));
    dc1_dD = (2*pi)^2*real(sum( ((-N/2:N/2)').^2 .* (conj(v) .* u - 2*conj(p).*T)));
    
    dc2_dF = -2*pi*imag(sum( (-N/2:N/2)' .* conj(T1) .*p));
    dc2_dD = (2*pi)^2 *real(sum( ((-N/2:N/2)').^2 .* conj(T1) .*p));
    
    dc1_dbr = zeros(size(u));
    dc1_dbi = zeros(size(u));
    dc2_dbr = zeros(size(u));
    dc2_dbi = zeros(size(u));
    offset = N/2 + 1;
   
    for idx = 1:length(dc1_dbr)
        k = idx - offset;
        l = -N/2:N/2;
        dxs = l((l - k >= -N/2) & (l-k <= N/2)) + offset;
        dc1_dbr(idx) = dc1_dbr(idx) + real(sum(conj(v(dxs)).* u(dxs -k) - 2* conj(p(dxs-k)).*T(dxs)));
        dc1_dbi(idx) = dc1_dbi(idx) - imag(sum(conj(v(dxs)).* u(dxs -k) + 2* conj(p(dxs-k)).*T(dxs)));
        
        dc2_dbr(idx) = dc2_dbr(idx) + real(sum(conj(T1(dxs)).*p(dxs-k)));
        dc2_dbi(idx) = dc2_dbi(idx) - imag(sum(conj(T1(dxs)).*p(dxs-k)));
        
        dxs = l((l + k >= -N/2) & (l+k <= N/2)) + offset;
        dc1_dbr(idx) = dc1_dbr(idx) + real(sum(conj(v(dxs)).* u(dxs +k) - 2* conj(p(dxs+k)).*T(dxs)));
        dc1_dbi(idx) = dc1_dbi(idx) + imag(sum(conj(v(dxs)).* u(dxs +k) + 2* conj(p(dxs+k)).*T(dxs)));
        
        dc2_dbr(idx) = dc2_dbr(idx) + real(sum(conj(T1(dxs)).*p(dxs+k)));
        dc2_dbi(idx) = dc2_dbi(idx) + imag(sum(conj(T1(dxs)).*p(dxs+k)));
    end
    
    c1 = T2 - t2; 
    dc1 =repack(M,dc1_dar + 1j*dc1_dai,reduce(dc1_dbr) + 1j *reduce(dc1_dbi) ,dc1_dD,dc1_dF);
    
    c2 = real(p(offset)) - 1/2;
    dc2 =repack(M,reduce(dc2_dar+1j*dc2_dai) ,reduce(dc2_dbr + 1j*dc2_dbi)  ,dc2_dD,dc2_dF);
    
    c=[];
    DC = [];
    ceq = [c1,c2];
    DCeq = [dc1,dc2];
    
end