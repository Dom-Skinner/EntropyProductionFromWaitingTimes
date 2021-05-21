function [a,c,Dk,Fk] = unpack(x,N,M)
    a_r = x(1:M);
    a_i = x(M+1:2*M);
    c_r = x(2*M+1:3*M+1);
    c_i = x(3*M+2:4*M+1);
    D_ = x(4*M+2);
    F_ = x(4*M+3);

    c = zeros(N/2+1,1);
    c(1) = c_r(1);
    c(2:M+1) = c_r(2:end) + 1j*c_i;

    a = zeros(N/2+1,1);
    a(1) = 1;
    a(2:M+1) = a_r + 1j*a_i;

    Dk = zeros(N/2+1,1);
    Dk(1) = D_;

    Fk = zeros(N/2+1,1);
    Fk(1) = F_;
end