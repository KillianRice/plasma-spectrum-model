function [tsol,psol] = ode_rk_integrator(t,p0,dpdt)
% t (vector double): time window (start and end) to return solutions for with units gam422^{-1}
% p0 (vector double): p0(k) is the ion state population of state |k> at the initial time, t(1).
% dpdt (afun @(t,p)): rate equations as an anonymous function of t and p (note information like
%                       scattering rate, decay rate, etc has already been fixed). This function is
%                       solely for solving the rate equations after they are defined.
% tsol (vector double): time points solutions are returned for with units gam422^{-1}
% psol (matrix double): p(k,n) is the population of state |k> at time tsol(n) returned by the differential
                % equation solver

h = 0.25; % initial timestep for ODE solver, handles initial rapid changes in population until quasi-equilibrium is reached.
maxiter = ceil(t(end)/h+1);

% preallocate space to tsol and psol for speed
tsol = zeros(1,maxiter);
psol = zeros(length(p0),ceil(t(end)/h));

% place initial conditions in solution arrays
tsol(1) = t(1);
psol(:,1) = p0;

% iteratively solve rate equations using Cash-Karp Runge-Kutta method.
iter = 1;
while tsol(iter) < t(end)
    iter = iter + 1;

    if tsol(iter-1) > 4 % increase timestep once initial rapid changes in population have subsided
        h = 2; % this is timestep for capturing optical pumping between ground states
    end

    [tsol(iter),psol(:,iter)] = rk38_method(tsol(iter-1),psol(:,iter-1),h,dpdt);
    
%     basic error checking
    eps = 1e-3;
    if max(psol(:,iter) < -eps)
        error('ODE solver failed.')
    end
end

% remove unused portion of preallocated solution arrays
tsol(iter+1:end) = [];
psol(:,iter+1:end) = [];

end
