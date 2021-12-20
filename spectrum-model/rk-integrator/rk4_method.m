function [t_np1,y_np1] = rk4_method(t_n,y_n,h,f)
% t_n (1x1 double): time at initial step
% t_np1 (1x1 double): time after rk4 step
% h (1x1 double): time step -> h = t_np1 - t_n
% y_n (1xm double): quantity of interest at initial time step
% y_np1 (1xm double): quantity of interest after initial time step, same size as y_np1
% f (anonymous fun @(t,y)): rate of change of y (i.e., f = dy/dt)

k1 = f(t_n,y_n);
k2 = f(t_n+h/2,y_n+h.*k1./2);
k3 = f(t_n+h/2,y_n+h.*k2./2);
k4 = f(t_n+h,y_n+h.*k3);

t_np1 = t_n + h;
y_np1 = y_n + h/6.*(k1+2.*k2+2.*k3+k4);

end
