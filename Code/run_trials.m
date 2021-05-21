function sig_est = run_trials(ntrials,sig_est,nvars,sig_inx,hess,rand_init,gs)
% Run a number of trials and update the matrix sig_est accordingly
ctol = 1e-7;
otol = 2e-7;
fvals = 2*ones(ntrials,1); % min is always <= 2
for j = 1:ntrials
    
    if ~rand_init
        wpvals = linspace(1.01,10,ntrials);
        wp = wpvals(j);
        B0 = -ones(nvars,nvars);
        for k= 1:size(B0,1)
            B0(k,k) = nvars+wp-1;
            if k < nvars
                B0(k,k+1) = -wp;
            end
        end
        B0 = B0/sum(B0(:));
        x0 = ones(nvars^2 + nvars,1);
        x0(1:nvars^2) = reshape(B0,nvars^2,1);
        x0(nvars^2 + 1 : end) = 1/nvars;
    else
         x0 = false;
    end
            
    
    [~,fval,exitflag] = min_via_fmin(nvars, sig_est(sig_inx,1),ctol,otol,x0,hess,gs);
    if exitflag >= 1 
        fvals(j) = fval;
    else
        disp('did not converge')
    end
end

if (min(fvals) <  sig_est(sig_inx,nvars)) || (sig_est(sig_inx,nvars) == 0)
    sig_est(sig_inx,nvars)  = min(fvals);
end

end
