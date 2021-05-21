function vq = gamma_precomp(x_in,invert)
% use invert = true if inputing in <t^2>/<t>^2 else use false if inputting
% sigma
load('nvars_comb.mat','sig_keep')
sig = sig_keep(:,1);
f = min(sig_keep(:,2:end),[],2);

if invert
    x_in(x_in > 2) = 2;
    vq = interp1(f,sig,x_in,'pchip');
else
    vq = interp1(sig,f,x_in,'pchip');
end
end