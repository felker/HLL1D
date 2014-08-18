function [rho_new,m_new,E_new] = riemann_hll(rho,u,p,m,E,gamma,Nx,lambda)

%Fluxes in conservation equation
f_rho = m;
f_u = rho.*u.*u + p;
f_p = u.*(1/2*rho.*u.*u + gamma/(gamma - 1)*p);
%these arent appropriate names since 
%since the fluxes are for m,E, not the primitives, u,p
%f_p = (m/rho)*(e+p) = u.*(rho.*eps + 1/2.*rho.*u^2 + p)
%Intercell fluxes
fc_rho = zeros(1,Nx);
fc_u = zeros(1,Nx);
fc_p = zeros(1,Nx);

for i=1:Nx-1
    [sl, sr] = wavespeed_estimator(rho,u,p,gamma,i,1);
    if sl>=0 
        fc_rho(1,i) = f_rho(1,i); 
        fc_u(1,i) = f_u(1,i); 
        fc_p(1,i) = f_p(1,i); 
    elseif sl<=0 && sr>=0 
        fc_rho(1,i) = (sr*f_rho(1,i) - sl*f_rho(1,i+1) + sl*sr*(rho(1,i+1)-rho(1,i)))/(sr-sl);
        fc_u(1,i) = (sr*f_u(1,i) - sl*f_u(1,i+1) + sl*sr*(m(1,i+1)-m(1,i)))/(sr-sl);
        fc_p(1,i) = (sr*f_p(1,i) - sl*f_p(1,i+1) + sl*sr*(E(1,i+1)-E(1,i)))/(sr-sl);
    elseif sr<=0
        fc_rho(1,i) = f_rho(1,i+1); 
        fc_u(1,i) = f_u(1,i+1); 
        fc_p(1,i) = f_p(1,i+1); 
    end
end
%update conserved variables
rho_new = rho + lambda*(circshift(fc_rho,[0,+1]) - fc_rho); %check bndry conditions
m_new = m + lambda*(circshift(fc_u,[0,+1]) - fc_u); %check bndry conditions
E_new = E + lambda*(circshift(fc_p,[0,+1]) - fc_p); %check bndry conditions
% i think circshift is fine because we throw out 1,1 anyway
end 

