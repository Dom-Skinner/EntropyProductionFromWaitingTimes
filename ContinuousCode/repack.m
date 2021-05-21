function x = repack(M,a,c,Dk,Fk)
    x = zeros(4*M+3,1);
    
    c_r = real(c(1:M+1));
    c_r(1) = 0.5*c_r(1);
    c_i = imag(c(2:M+1));
    
    a_r = real(a(2:M+1));
    a_i= imag(a(2:M+1));

    x(1:M) = a_r;
    x(M+1:2*M) =  a_i;
    x(2*M+1:3*M+1) = c_r;
    x(3*M+2:4*M+1) = c_i ;
    x(4*M+2) = Dk(1);
    x(4*M+3) =  Fk(1);

end