function [A,C] = jacobians(xuno,xdos,current,gamma, alfa, v0, vl, beta, dt, Ecrit)
    dv_dx1 = -current;
    if xdos==0
        dv_dx2 = 0;
    else 
        dv_dx2 = gamma*(v0-vl)*exp(gamma*(xdos-1)) + alfa*vl + ...
            (beta*(1-alfa)*vl/(2*sqrt(xdos)))*exp(-beta*sqrt(xdos));
    end
    factor = current*dt/Ecrit;
    
    A = [1,0; factor*dv_dx1, 1-factor*dv_dx2];
    C = [dv_dx1, dv_dx2];
end