%1D Sod Shock tube test
clear all;
Nx = 100;
Lx = 1.;
dx = Lx/Nx;
x = linspace(0.0,1.0,Nx-2);

nsteps = 50;

%initial condition: (physical variables are 2:Nx-1)
% primitive variables, left side of tube
rho(1,1:Nx/2) = 1.0;
u(1,1:Nx) = 0.0;
p(1,1:Nx/2) = 1.0;

%primitive variables, right side of tube
rho(1,Nx/2+1:Nx) = 0.125;
p(1,Nx/2+1:Nx) = 0.1;

%boundary conditions?
%conserved quantities in time derivative:
%rho
%rho*u = m
%e= rho(0.5*u*u + eps), closed with EOS p = (gamma -1)*rho*eps -->
%e = p/(gamma-1) +0.5*u*u*rho
gamma = 1.4; 
m = rho.*u;
E = p/(gamma-1) + 0.5*u.*u.*rho;

%CFL parameter
cfl =0.8;

t=0.0;
%dt=0.0; 
%dt =0.004111 %this is for Davis88 reproduction
for i=1:nsteps
if mod(i,10) ==0
    plot(x,u(2:Nx-1),'-o'); 
end
%REA Method
%Reconstruction: Use piecewise constant, first order reconstruction
dt = calculate_dt(rho,u,p,gamma,dx,cfl);
  
%Solve Riemann Problems
[rho_new,m_new,E_new] = riemann_hll(rho,u,p,m,E,gamma,Nx,dt/dx);

%apply transmission boundary conditions
rho_new(1,1) = rho_new(1,2);
m_new(1,1) = m_new(1,2);
E_new(1,1) = E_new(1,2);
rho_new(1,Nx) = rho_new(1,Nx-1);
m_new(1,Nx) = m_new(1,Nx-1);
E_new(1,Nx) = E_new(1,Nx-1);

rho = rho_new;
m = m_new;
E = E_new;
%update primitive variables from conserved variables
u = m./rho;
p = (gamma-1).*(E - 0.5*m.*m./rho); 

t= t+dt
end
