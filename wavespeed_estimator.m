%As in Davis88 and Einfeldt, this completes the approxmiation for Harten,
%Lax and van Leer's approxmiate Riemann solvers
function [sl,sr] = wavespeed_estimator(rho,u,p,gamma,i,METHOD)
rho_l = rho(i);
rho_r = rho(i+1);
ul = u(i);
ur = u(i+1);
pl = p(i);
pr = p(i+1);

%enthalpy derivation
%recall, closure is ideal EOS p = (gamma - 1)*eps*rho
%eps = total energy per unit mass
%e = total energy per unit length = rho*eps + 0.5*rho*u^2
%H = total enthalpy, rho*H = e + p = rho*eps +1/2*rho*u^2 + p
%H = eps +0.5*u^2 + p/rho = p/((gamma-1)*rho) +0.5*u^2 + p/rho 
H_l = gamma*pl/((gamma-1)*rho_l) + 0.5*ul*ul;
H_r = gamma*pr/((gamma-1)*rho_r) + 0.5*ur*ur;
%H_l = pl/((gamma-1)*rho_l) + 0.5*ul*ul + pl/rho_l;

switch METHOD
    case 1 %Davis, unbounded eqn 4.6 for HLL
        %Speed of sound in the left cell
        c_l = sqrt(gamma*pl/rho_l);
        c_r = sqrt(gamma*pr/rho_r); 
        %these are the first and last eigenvalues of A(u)
        sl = ul - c_l;  
        sr = ur + c_r;
        %Using the Roe linearized equaiton, Einfeldt88 has a different wavespeed
%estimator. Davis acknowledges its existence but doesnt implement it. 
    case 2 %Einfeldt Roe-linearized for HLLE, eqn 5.1 in Einfeldt88
        %averaged velocity 
        u_half = (sqrt(rho_l)*ul + sqrt(rho_r)*ur)/(sqrt(rho_l)+sqrt(rho_r));
        %averaged total enthalpy
        H_half = (sqrt(rho_l)*H_l + sqrt(rho_r)*H_r)/(sqrt(rho_l) + sqrt(rho_r));
        %averaged speed of sound
        c_half = sqrt((gamma-1)*(H_half - 1/2*u_half*u_half));
        sl = u_half - c_half;
        sr = u_half + c_half;
end
%An alternative would be to use pressure based estimates-- See Toro

end