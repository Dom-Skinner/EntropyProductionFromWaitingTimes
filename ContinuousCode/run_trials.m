function sig = run_trials(N,M,ntrials,t2)

otol = 1e-6;
ctol = 1e-6;
sig = -2*ones(ntrials,1);

for i = 1:ntrials
    [~,~,~,~,~,~,fval,exitflg] = min_cont(N,M, t2,otol,ctol);
    if exitflg > 0
        sig(i) = fval;
    end
end

sig = min(sig(sig>0));
end