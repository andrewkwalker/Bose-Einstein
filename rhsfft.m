function rhs_vec = rhsfft(t,psi_vec_f,K,A1,A2,A3,B1,B2,B3,X,Y,Z,n)

    psi_f = reshape(psi_vec_f,n,n,n);
    psi = ifftn(psi_f);
    
    rhs = 1i*(.5*-K.*psi_f - ...
              fftn((abs(psi).^2).*psi) + ...
              fftn((A1*((sin(X)).^2)+B1).*...
                   (A2*((sin(Y)).^2)+B2).*...
                   (A3*((sin(Z)).^2)+B3).*psi));
          
    rhs_vec = reshape(rhs,n^3,1);

end