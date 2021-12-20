function [t_np1,y_np1] = rk38_method(t_n,y_n,h,f)
% t_n (1x1 double): time at initial step
% t_np1 (1x1 double): time after rk4 step
% h (1x1 double): time step -> h = t_np1 - t_n
% y_n (1xm double): quantity of interest at initial time step
% y_np1 (1xm double): quantity of interest after initial time step, same size as y_np1
% f (anonymous fun @(t,y)): rate of change of y (i.e., f = dy/dt)
a2 = [1/3];
a3 = [-1/3 1];
a4 = [1 -1 1];
b = [1/8 3/8 3/8 1/8];
c = [0 1/3 2/3 1];

k1 = f(t_n,y_n);
k2 = f(t_n+h*c(2),y_n+h.*(a2(1).*k1));
k3 = f(t_n+h*c(3),y_n+h.*(a3(1).*k1+a3(2).*k2));
k4 = f(t_n+h*c(4),y_n+h.*(a4(1).*k1+a4(2).*k2+a4(3).*k3));

t_np1 = t_n + h;
y_np1 = y_n + h.*(b(1).*k1+b(2).*k2+b(3).*k3+b(4).*k4);

end
