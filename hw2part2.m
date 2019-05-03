
alphaD = 1; % alpha_d
deltaJMin = 1; % min change in cost
currTime = 0; % t_curr
defaultControlDuration = 0.1; % deltaT_init
nominalControlU1 = 0; % nominal control u1
omega = 0.5; % scale factor omega
predictionHorizon = 0.2; % T
samplingTime = 1; % t_s
tCalcMax = 1; % t_calc max time for iterative control calculations
maxBacktrackingIters = 2; % k_max

%%% want to find a control vector in each time window
%%%  then find the time tau at which to apply that control vector
%%%  then find the application duration lambda

while currTime < inf
    currTime = 0
    
end


%%%
function j = J(vals)
  % vals in form of:
  % [[x10, x11, ...., x1N-1],
  %  [x20, x21, ...., x2N-1],
  %  [x30, x31, ...., x3N-1],
  %  [u10, u11, ...., u1N-1],
  %  [u20, u21, ...., u2N-1]]
  global N Q R P1
  xs = vals(1:3, 1:N);
  us = vals(4:5, 1:N);

  sumsofar = 0;
  for i = 1:N-1
    xi = [xs(1, i); xs(2, i); xs(3, i)];
    ui = [us(1, i); us(2, i)];
    xdi = Fdest(i);
    xidiff = xi - xdi;
    sumsofar = sumsofar + transpose(xidiff) * Q * xidiff + transpose(ui) * R * ui;
  end
  sumsofar;
  xsn = [xs(1, N); xs(2, N); xs(3, N)];
  xdn = Fdest(N);
  xndiff = xsn - xdn;
  sumsofar = 0.5*sumsofar + 0.5*(transpose(xndiff) * P1 * xndiff);
  j = sumsofar;
end
function dj = DJ(xs, controls, zs, vs)

  global N Q R P1

  sumsofar = 0;
  for i = 1:N-1
    xi = [xs(1, i); xs(2, i); xs(3, i)];
    ci = [controls(1, i); controls(2, i)];
    zsi = zs(i,:);
    vsi = vs(:,i);
    xdi = Fdest(i);
    xidiff = xi - xdi;
    sumsofar = sumsofar + transpose(xidiff) * Q * transpose(zsi) + transpose(ci) * R * vsi;
  end
  sumsofar;
  xsn = [xs(1, N); xs(2, N); xs(3, N)];
  xdn = Fdest(N);
  xndiff = xsn - xdn;
  sumsofar = 0.5*sumsofar + 0.5*(transpose(xndiff) * P1 * transpose(zs(N,:)));
  dj = sumsofar;
end
%%%
function xdest = Fdest(i)
  global dtval;

  xdest = [dtval * 2 * i / pi; 0; pi/2];
end