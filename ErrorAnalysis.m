clf;
clc;
clear All;
fig = 3;
format long
dx=.12
%% Heat Equation
for n=0:8
  load(sprintf('Heat_NoFlux_dt%d.mat', n));
  %times(n+1) = log2(dt);
  times(n+1) = dt;
  
  %% Makes sense to average the norms by number of points
  mimeticError(n+1) = (norm(exact-mimetic, 2));
  finiteDifferenceError(n+1) = (norm(exact-finiteDifference, 2));
  
%  mimeticError(n+1)= dx.^2*sum(abs(exact(:)-mimetic(:)).^1);
%   
%  finiteDifferenceError(n+1)=dx.^2*sum(abs(exact(:)-finiteDifference(:)).^1);
  
  %times=fliplr(times);
  
    
end
 
 mimeticError=fliplr(mimeticError)
 finiteDifferenceError=fliplr(finiteDifferenceError)

figure(fig); fig = fig + 1;
%hold on;
loglog(times, finiteDifferenceError, 'k',times, mimeticError, 'ro--');
xlabel('dt'); ylabel('Error');
title('Error Comparison for Heat Eq. w/ No Flux B.Cs');
legend('FD','Mimetic');
times