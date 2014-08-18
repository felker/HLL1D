function dt = calculate_dt(rho,u,p,gamma,dx,cfl)
  % value of local flow speed + local sound speed
  u_local = abs(u) + sqrt(gamma*p./rho);
  u_max = max(u_local)
  max(abs(u))
  dt = cfl*dx/u_max;
end