function x = s_extrap(F,X)
% Function to extrapolate to n = \infty
ObjectiveFunction = @(x) fit_res(x,F,X);
C1 = linspace(1,min(F),1000);
[R,~] = fit_res(C1,F,X);
idx = find( (R(1:end-2) > R(2:end-1)) & (R(3:end) > R(2:end-1)));
if isempty(idx)
    x = min(F);
    return 
end
lb = C1(idx);
ub = C1(idx+2);


options = optimoptions('fmincon', 'OptimalityTolerance',1e-6,...
    'MaxFunctionEvaluations', 1.2e+02,'MaxIterations',3e2,'SpecifyObjectiveGradient',true);
options.Display = 'iter';
[x,~,~,~]  = fmincon(ObjectiveFunction,0.5*(lb+ub),[],[],[],[],lb,ub,[],options);



function [R,dR] = fit_res(C1,F,X)
    if length(C1) > 1
        R = zeros(size(C1));
        dR = zeros(size(C1));
        for i = 1:length(C1)
            [R(i),dR(i)] = fit_res(C1(i),F,X);
        end
    else
        Y = log(F-C1);
        dY = - 1./(F - C1);
        Y = Y - mean(Y);
        dY = dY - mean(dY);
        X = X - mean(X);

        C3 = sum(X.*Y)/sum(X.^2);
        dC3 = sum(X.*dY)/sum(X.^2);

        R = sum( (Y - C3*X).^2);
        dR = sum( 2*(dY - dC3*X).*(Y - C3*X));
    end
end

end

