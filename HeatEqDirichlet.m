N = 50; % number of grid points

%% X Domain
a = 0;
b = 2*pi;

%% X discretization
dx = (b-a)/N;
dx2 = dx/2;
x = [a (a+dx2):dx:(b-dx2) b]';

%% Time Domain
t0 = 0;
tf = 5;
dt = (dx^2)/(4*alpha);
t = [t0:dt:tf];

%% Matlab Solver Approach
sol = pdepe(0, @pdefun, @icfun, @bcfun, x, t);

figure(8);
plot(x, sol(end,:), 'o-');
str = sprintf('Matlab');
title(str)
xlabel('x')
ylabel('dff')
xlim([a b]);
ylim([-2 2]); 

figure(9);
mesh(x, t, sol);
title("matlab");
xlabel('x'); ylabel('t'); zlabel('u');

function [c,f,s] = pdefun(x, t, u, DuDx)
  coeff= (1/.15).^2;  %Coefficient needed for Allen-Cahn PDE.
  c = 1;
  f = DuDx;
  s = -coeff*(u.^3-u);
end

function u0 = icfun(x)
  u0= 1.2*(rand(size(x))-1*rand(size(x)));
  %% u0 = 0.5*(1+sin(2*2*3.14 * x/(2*pi)));
end

function [pl,ql,pr,qr] = bcfun(xl, ul, xr, ur, t)
  pl = ul;
  ql = 0;
  pr = ur;
  qr = 0;
end

  

