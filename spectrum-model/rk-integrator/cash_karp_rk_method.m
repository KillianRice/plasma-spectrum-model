function [t_np1,y_np1] = cash_karp_rk_method(t_n,y_n,h,f)
% t_n (double): time at initial step
% t_np1 (double): time after rk4 step
% h (double): time step -> h = t_np1 - t_n
% y_n (vector double): quantity of interest at initial time step
% y_np1 (vector double): quantity of interest after initial time step, same size as y_np1
% f (anonymous fun @(t,y)): rate of change of y (i.e., f = dy/dt)

a2 = 1/5;
a3 = [3/40 9/40];
a4 = [3/10 -9/10 6/5];
a5 = [-11/54 5/2 -70/27 35/27];
a6 = [1631/55296 175/512 575/13824 44275/110592 253/4096];
b = [2825/27648 0 18575/48384 13525/55296 277/14336 1/4];
c = [0 1/5 3/10 3/5 1 7/8];

k1 = f(t_n,y_n);
k2 = f(t_n+h*c(2),y_n+h.*(a2(1).*k1));
k3 = f(t_n+h*c(3),y_n+h.*(a3(1).*k1+a3(2).*k2));
k4 = f(t_n+h*c(4),y_n+h.*(a4(1).*k1+a4(2).*k2+a4(3).*k3));
k5 = f(t_n+h*c(5),y_n+h.*(a5(1).*k1+a5(2).*k2+a5(3).*k3+a5(4).*k4));
k6 = f(t_n+h*c(6),y_n+h.*(a6(1).*k1+a6(2).*k2+a6(3).*k3+a6(4).*k4+a6(5).*k5));

t_np1 = t_n + h;
y_np1 = y_n + h.*(b(1).*k1+b(2).*k2+b(3).*k3+b(4).*k4+b(5).*k5+b(6).*k6);

end
