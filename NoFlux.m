%% Allen-Cahn and Heat Equations with No-Flux Boundary conditions.
addpath('./mole_MATLAB/');
run('Configuration.m');

%% Domain
X = [a (a+dx2):dx:(b-dx2) b]';
T = [t0:dt:(ceil(tf/dt)*dt)];

%% Initial Condition
ICfun = @(X) 1.2*(rand(size(X))-1*rand(size(X)));
ICfun = @(X) 0.15*(1.0+sin(3*2*3.14/(b-a) * X));
ICfun = @(X) cos(X);
u0 = ICfun(X);

%% "Standard" Finite Differences Approach
finiteDifference = zeros(length(T), length(X)); % to store solutions

%% Time discretization
t = t0;
timesteps = 1;

%% Initial Condition
u = u0;
unew = u0;
finiteDifference(1,:) = u0;

%% Finite Difference Operator Matrix
%% N
v = dt/dx^2;
a0 = v*ones(N+2, 1);
a1 = -2*v*ones(N+2, 1);
A = spdiags([a0 a1 a0], [-1 0 1], N+2, N+2);

%% Use Saulo's non-uniform spacing scheme near edges:
%% X = a + dx/2 
A(2, 1)   = (8/3)*v; A(2, 2)     = -4*v; A(2, 3)   = (4/3)*v;
%% X = b - dx/2 
A(N+1, N) = (4/3)*v; A(N+1, N+1) = -4*v; A(N+1, N+2) = (8/3)*v;

%% Enforce No-Flux B.Cs:
%% Really we are using the same centered difference laplacian over the boundary,
%% however we treat the A(1,0) point as equal to A(1,2), and likewise
%% for the other edge.
vedge = 4*v;
%% X = a
A(1, 1) = -2*vedge; A(1, 2) = 2*vedge;
%% X = b
A(N+2, N+1) = 2*vedge; A(N+2, N+2) = -2*vedge;

if (explicit)
  FD = speye(size(A)) + alpha*A;
else
  FD = speye(size(A)) - alpha*A;
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

    timesteps = timesteps+1;
    finiteDifference(timesteps,:) = u;
end

u_finiteDifference = u;

%% Mimetic Finite Differences Approach
mimetic = zeros(length(T), length(X));

%% Time discretization
t = t0;
timesteps = 1;

%% Initial Condition
u = u0;
unew = u0;
mimetic(1,:) = u0;

%% Mimetic Operator Matrix
order = 2;

L = lap(order, N, dx);
R = robinBC(order, N, dx, 1, 1);

if (explicit)
  MFD = speye(size(L)) + alpha*dt*(L-R);
else
  MFD = speye(size(L)) - alpha*dt*(L-R);
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

    u = unew;

    timesteps = timesteps+1;
    mimetic(timesteps,:) = u;
end

u_mimetic = u;

%% Matlab Solver
matlab = pdepe(0, @(x, t, u, DuDx) pde(x, t, u, DuDx, equation, loss_coeff, allen_coeff), @(x) u0(find(X==x)), @bcfun, X, T);

%% Exact Solutions
exact = zeros(length(T), length(X));
has_exact = false;
if (equation == HEAT)
  has_exact = true;
  for t=1:length(T)
    exact(t, :) = (exp(-T(t))*cos(X))';
  end
elseif (equation == HEATLL)
  has_exact = true;
  for t=1:length(T)
    exact(t, :) = exp(-T(t)*(1+loss_coeff))*cos(X)';
  end
end

%% Otherwise, use a reference
reference = zeros(length(T), length(X)); % to store solutions
has_exact = false;
if (do_plots && ~has_exact)
  t = 0;
  timesteps = 1;
  
  dt = dt / tSubd;

  unew = u0;
  reference(1,:) = u0;

  dx = dx / xSubd;
  Xsub = [a:(dx):b]';
  Nsub = length(Xsub);
  Indices = [1 (1+xSubd/2):(xSubd):(Nsub-xSubd/2) Nsub];
  u = ICfun(Xsub);
  unew = u;
  reference(1, :) = u(Indices);

  v = dt/dx^2;
  a0 = ones(Nsub, 1);
  a1 = -2*ones(Nsub, 1);
  A = spdiags([a0 a1 a0], [-1 0 1], Nsub, Nsub);
  A(1, 2) = 2;
  A(end, end-1) = 2

  if (explicit)
    RFD = speye(size(A)) + alpha*v*A;
  else
    RFD = speye(size(A)) - alpha*v*A;
  end

  while (timesteps < length(T))
    for i=1:tSubd
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
    reference(timesteps, :) = u(Indices);
  end
end

%% Plots
if (do_plots)
  figure(figures_so_far); figures_so_far = figures_so_far + 1;
  mesh(X,T,finiteDifference);
  title("Finite Difference");
  xlabel('x'); ylabel('t'); zlabel('u');

  figure(figures_so_far); figures_so_far = figures_so_far + 1;
  mesh(X,T,mimetic);
  title("Mimetic");
  xlabel('x'); ylabel('t'); zlabel('u');

  figure(figures_so_far); figures_so_far = figures_so_far + 1;
  mesh(X,T,matlab);
  title("Matlab");
  xlabel('x'); ylabel('t'); zlabel('u');

  if (has_exact)
    compareTo = exact;
    compareLabel = "Exact";
  else
    compareTo = reference;
    compareLabel = "Reference";
  end
  
  figure(figures_so_far); figures_so_far = figures_so_far + 1;
  mesh(X,T,abs(compareTo-finiteDifference));
  title(sprintf("Finite Difference vs %s", compareLabel));
  xlabel('x'); ylabel('t'); zlabel('diff');

  figure(figures_so_far); figures_so_far = figures_so_far + 1;
  mesh(X,T,abs(compareTo-matlab));
  title(sprintf("Matlab vs %s", compareLabel));
  xlabel('x'); ylabel('t'); zlabel('diff');

  figure(figures_so_far); figures_so_far = figures_so_far + 1;
  mesh(X,T,abs(compareTo-mimetic));
  title(sprintf("Mimetic vs %s", compareLabel));
  xlabel('x'); ylabel('t'); zlabel('diff');

  figure(figures_so_far); figures_so_far = figures_so_far + 1;
  mesh(X,T,abs(finiteDifference-mimetic));
  title("Finite Difference vs Mimetic");
  xlabel('x'); ylabel('t'); zlabel('diff.');
  
  figure(figures_so_far); figures_so_far = figures_so_far + 1;
  clf; hold on;
  plot(X,u0, "k:", 'DisplayName', 'I.C');
  plot(X, matlab(end,:), "k:", 'DisplayName', 'Matlab');
  plot(X, mimetic(end,:), "r-", 'DisplayName', 'Mimetic');
  plot(X, finiteDifference(end,:), "b--", 'DisplayName', 'Finite Difference');
  plot(X, compareTo(end,:), "k.", 'DisplayName', compareLabel);
  title(sprintf("Final State (%s, dt=%g)", eqTitle, tf));
  xlabel('x');
  ylabel('u');
  xlim([a b]);
  ylim([-2 2]);
  legend();
  hold off;

  figure(figures_so_far); figures_so_far = figures_so_far + 1;
  clf; hold on;
  %% plot(X, abs(matlab(end,:)-compareTo(end,:)), "k:", 'DisplayName', 'Matlab');
  plot(X, abs(mimetic(end,:)-compareTo(end,:)), "r-", 'DisplayName', 'Mimetic');
  plot(X, abs(finiteDifference(end,:)-compareTo(end,:)), "b--", 'DisplayName', 'Finite Difference');
  title(sprintf("Final Error (%s, dt=%g)", eqTitle, tf));
  xlabel('x');
  ylabel('u');
  xlim([a b]);
  legend();
  hold off;
end

%% Save Outputs
if (save_outputs)
  save(sprintf('./data/%s_NoFlux_dt%d_dx%d.mat', eqLabel, dtMagnitude, dxMagnitude), "dx", "dt", "exact", "mimetic", "finiteDifference", "matlab", "reference");
end

%% Function for matlab solver
function [c,f,s] = pde(x, t, u, DuDx, equation, loss_coeff, allen_coeff)
  c = 1;
  f = DuDx;
  switch (equation)
    case 0
      s = 0;
    case 1
      s = -loss_coeff*u;
    case 2
	  s = allen_coeff*(u-u.^3);
  end        
end

function [pl,ql,pr,qr] = bcfun(xl, ul, xr, ur, t)
  pl = 0;
  ql = 1;
  pr = 0;
  qr = 1;
end

    
