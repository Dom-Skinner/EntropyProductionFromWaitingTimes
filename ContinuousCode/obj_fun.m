function [f,df] = obj_fun(x,N,M)

    [ak,ck,Dk,Fk] = unpack(x,N,M);
    L = L_matrix(Fk,Dk,ck);
    p = -0.5*reduce((L') \ make_full(ak));
    
    f = entropy_rate(ak,ck,p,Fk,Dk);
    [P,Q,R] = entropy_grad(ak,ck,p);
    p = make_full(p);


    w = conj(L \ conj(P));
    ds_dF  = real(Fk(1)/Dk(1) + 2*pi* 1j * sum( (-N/2 : N/2)' .* w.*p));
    ds_dD  = real(-0.5*Fk(1)^2/Dk(1)^2 + 4*pi^2 * sum( ((-N/2 : N/2).^2)' .* w .*p));
    ds_dar = real(-w + Q);
    ds_dai = imag(w - Q);
    ds_dbr = real(R);
    ds_dbi = -imag(R);
    offset = N/2 + 1;
  
    for idx = 1:length(ds_dbi)
        k = idx - offset;
        l = -N/2:N/2;
        dxs = l((l - k >= -N/2) & (l-k <= N/2)) + offset;
        ds_dbr(idx) = ds_dbr(idx) + real(sum(w(dxs).* p(dxs -k)));
        ds_dbi(idx) = ds_dbi(idx) - imag(sum(w(dxs).*p(dxs -k)));
    
        dxs = l((l + k >= -N/2) & (l+k <= N/2)) + offset;
        ds_dbr(idx) = ds_dbr(idx) + real(sum(w(dxs).* p(dxs +k)));
        ds_dbi(idx) = ds_dbi(idx) + imag(sum(w(dxs).*p(dxs +k)));
       
    end
    
    ds_dar = reduce(ds_dar);
    ds_dai = reduce(ds_dai);
    ds_dbr = reduce(ds_dbr);
    ds_dbi = reduce(ds_dbi);
    
    df = repack(M,ds_dar + 1j*ds_dai,ds_dbr + 1j*ds_dbi,ds_dD,ds_dF);

    
function f = make_full(fk)
    f = [conj(fk(end:-1:1)); fk(2:end)];
end

function f = reduce(fk)
    f = fk((numel(fk)-1)/2+1 :end);
end

end