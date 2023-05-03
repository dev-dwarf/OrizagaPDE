%% Allen-Cahn and Heat Equations with Periodic Boundary conditions.
addpath('./mole_MATLAB/');
run('Configuration.m');

%% Domain
dx = dx / (2^dxMagnitude);
dx2 = dx / 2;
X = [a (a+dx2):dx:(b-dx2) b]';
Xsim = X(1:end-1);
Xdis = [X; (b-a)+X(2:end);];

dt = (dx2^2)/2;
dt = dt / (2^dtMagnitude);
T = [t0:dt:(floor(tf/dt)*dt)];

%% Initial Condition
u0 = 1.2*(rand(size(Xsim))-1*rand(size(Xsim)));
u0 = 0.5*(sin(2*3.14/(b-a) * Xsim));

%% "Standard" Finite Differences Approach
finiteDifference = zeros(length(T), length(Xdis)); % to store solutions

%% Time discretization
timesteps = 1;

%% Initial Condition
u = u0;
unew = u0;
U = [u0; u0; u0(1)];
finiteDifference(1,:) = U;

%% Finite Difference Operator Matrix
v = alpha*dt/dx^2;
a0 = v*ones(N+2, 1);
a1 = -2*v*ones(N+2, 1);
A = spdiags([a0 a1 a0], [-1 0 1], N+2, N+2);

%% Use Saulo's non-uniform spacing scheme at boundaries
%% X = a + dx/2
A(2, 1)   = (8/3)*v; A(2, 2)     = -4*v; A(2, 3)   = (4/3)*v;
%% X = b - dx/2
A(N+1, N) = (4/3)*v; A(N+1, N+1) = -4*v; A(N+1, N+2) = (8/3)*v;

%% Enforce Periodic B.Cs:
%% Use centered difference laplacian.
vedge = 4*v;
%% X = a
A(1, N+1) = vedge;       A(1, 1)     = -2*vedge; A(1, 2)   = vedge;
%% Apply the rule that U(a) == U(b), and drop last row/column.
A(:, 1) = A(:, 1) + A(:, end);
A = A(1:end-1,1:end-1);

if (explicit)
  FD = speye(size(A)) + A;
else
  FD = speye(size(A)) - A;
end

%% Integrate
while (timesteps < length(T))
    %% update U
    if (explicit)
      unew = FD * u;
    else
      unew = FD \ u;
    end

    switch (equation)
      case 1
        unew = unew - loss_coeff*dt*u;
      case 2
        unew = unew - allen_coeff*dt*(u.^3-u);
    end        

    u = unew;

    %% timestep
    timesteps = timesteps+1;
    U = [u; u; u(1)];
    finiteDifference(timesteps,:) = U;
end

u_finiteDifference = u;

%% Mimetic Finite Differences Approach
mimetic = zeros(length(T), length(Xdis));

%% Time discretization
timesteps = 1;

%% Initial Condition
u = u0;
unew = u0;
U = [u0; u0; u0(1)];
mimetic(1,:) = U;

%% Mimetic Operator Matrix
order = 2;

%% Instead of using MOLE's lap(), construct laplacian from div and grad ourselves.
%% This lets us impose boundary conditions on D and G as needed.
D = div(order, N, dx); % 1D Mimetic divergence operator
G = grad(order, N, dx);

%% Periodic BC imposed on the divergence operator
%% Taken from MOLE example code: examples_MATLAB/hyperbolic1D.m
%% https://github.com/jcorbino/mole/blob/686bfe038c5717957b67496dc11c552f872609f2/examples_MATLAB/hyperbolic1D.m
D(1,2) = 1/(2*dx);
D(1,end-1) = -1/(2*dx);
D(end,2) = 1/(2*dx);
D(end,end-1) = -1/(2*dx);

%% https://www.sciencedirect.com/science/article/pii/S0377042715005051#s000050

%% Construct Mimetic Laplacian
L = D*G;

%% Apply the rule that U(a) == U(b), and drop last row/column.
L(:, 1) = L(:, 1) + L(:, end);
L = L(1:end-1,1:end-1);

if (explicit)
  MFD = speye(size(L)) + alpha*dt*L;
else
  MFD = speye(size(L)) - alpha*dt*L;
end

%% Integrate
while (timesteps < length(T))
    %% update U
    if (explicit)
      unew = MFD * u;
    else
      unew = MFD \ u;
    end

    switch (equation)
      case 1
        unew = unew - loss_coeff*dt*u;
      case 2
        unew = unew - allen_coeff*dt*(u.^3-u);
    end        

    %% timestep
    timesteps = timesteps+1;
    u = unew;

    U = [u; u; u(1)];
    mimetic(timesteps,:) = U;
end

u_mimetic = u;

%% Reference Solution - Finite Difference, small timestep
reference = zeros(length(T), length(Xdis)); % to store solutions
timesteps = 1;
u = u0;
unew = u0;
U = [u0; u0; u0(1)];
reference(1,:) = U;

subdivide = 100;
dt = dt / subdivide;

if (explicit)
  RFD = speye(size(A)) + A / subdivide;
else
  RFD = speye(size(A)) - A / subdivide;
end

while (timesteps < length(T))
  for i=1:subdivide
    if (explicit)
      unew = RFD * u;
    else
      unew = RFD \ u;
    end

    switch (equation)
      case 1
        unew = unew - loss_coeff*dt*u;
      case 2
        unew = unew - allen_coeff*dt*(u.^3-u);
    end        

    u = unew;
  end
  
  timesteps = timesteps+1;
  U = [u; u; u(1)];
  reference(timesteps,:) = U;
end
  
%% Plots
if (do_plots)
  figure(figures_so_far); figures_so_far = figures_so_far + 1;
  mesh(Xdis,T,finiteDifference);
  title("Finite Difference");
  xlabel('x'); ylabel('t'); zlabel('u');

  figure(figures_so_far); figures_so_far = figures_so_far + 1;
  mesh(Xdis,T,mimetic);
  title("Mimetic");
  xlabel('x'); ylabel('t'); zlabel('u');

  figure(figures_so_far); figures_so_far = figures_so_far + 1;
  mesh(Xdis,T,reference);
  title("Reference");
  xlabel('x'); ylabel('t'); zlabel('u');

  figure(figures_so_far); figures_so_far = figures_so_far + 1;
  mesh(Xdis,T,abs(finiteDifference-reference));
  title("Finite Difference vs Reference");
  xlabel('x'); ylabel('t'); zlabel('diff.');

  figure(figures_so_far); figures_so_far = figures_so_far + 1;
  mesh(Xdis,T,abs(mimetic-reference));
  title("Mimetic vs Reference");
  xlabel('x'); ylabel('t'); zlabel('diff.');

  figure(figures_so_far); figures_so_far = figures_so_far + 1;
  clf; hold on;
  plot(Xdis, finiteDifference(1,:), "k:", 'DisplayName', 'I.C');
  plot(Xdis, mimetic(end,:), "r-", 'DisplayName', 'Mimetic');
  plot(Xdis, finiteDifference(end,:), "b--", 'DisplayName', 'FiniteDifference');
  plot(Xdis, reference(end,:), "k.", 'DisplayName', 'Reference');
  title(sprintf("Final State (%s, dt=%g)", eqTitle, tf));
  xlabel('x')
  ylabel('u')
  xlim([a b+b]);
  ylim([-2 2]);
  legend();
  hold off;

  figure(figures_so_far); figures_so_far = figures_so_far + 1;
  clf; hold on;
  plot(Xdis, abs(mimetic(end,:)-reference(end,:)), "r-", 'DisplayName', 'Mimetic');
  plot(Xdis, abs(finiteDifference(end,:)-reference(end,:)), "b--", 'DisplayName', 'FiniteDifference');
  title(sprintf("Final Error (%s, dt=%g)", eqTitle, tf));
  xlabel('x')
  ylabel('u')
  xlim([a b+b]);
  legend();
  hold off;
end

%% Save Outputs
if (save_outputs)
  save(sprintf('./data/%s_NoFlux_dt%d_dx%d.mat', eqLabel, dtMagnitude, dxMagnitude), "dx", "dt", "mimetic", "finiteDifference", "reference");
end
