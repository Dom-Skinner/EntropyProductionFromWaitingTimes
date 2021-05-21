function [x,B] = rand_init_cond(n_int,sig)

A = -rand(n_int+1,n_int+1);
A(1:n_int+2:end) = -A(1:n_int+2:end);
A(2:n_int+2:end) = 50*A(2:n_int+2:end);

f = @(A,s) 0.5*sum(sum( (A' - A) .* log( A./A'))) - s;
gf = @(A) A'./A - 1 - log(A./A');


A(end,end) = 1;

for i = 1:250
    % project column cond
    A = A - (A * ones(n_int+1,1)) * ones(1,n_int+1)/(n_int+1);
    % satisfy sign constraint
    A(A>0) = -A(A>0); 
    A(1:n_int+2:end) = -A(1:n_int+2:end);
    % project row cond
    A = A - ones(n_int+1,1)* (ones(1,n_int+1)*A)/(n_int+1);
    % satisfy sign constraint
    A(A>0) = -A(A>0);
    A(1:n_int+2:end) = -A(1:n_int+2:end);
    % timing constraint
    A(end,end) = 1;
    
    % Newton iterate for entropy
    lam = - f(A,sig)/ norm(gf(A))^2;
    A = A + lam*gf(A);
    A(A>0) = -A(A>0);
    A(1:n_int+2:end) = -A(1:n_int+2:end);
end


B = A(1:n_int,1:n_int);
x = rand(n_int,1);
x = x/sum(x);
if (min(sum(B)) < 0 ) || (min(sum(B')) < 0 )
    disp('recursing')
    [x,B] = rand_init_cond(n_int,sig/2); % doesn't converge - try with easier sig
end
end

