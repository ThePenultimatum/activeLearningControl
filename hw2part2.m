
clear
global N alphaD deltaJMin currTime defaultControlDuration nominalControlU1 omega predictionHorizon samplingTime tCalcMax maxBacktrackingIters endTime Q R P1 dtval totalN maxU1 maxU2 minU1 minU2 numItersControlDurationDefault

N = 100;

maxU1 = 5;
maxU2 = 5;
minU1 = 5;
minU2 = 5;

endTime = 2*pi;

alphaD = 1; % alpha_d desired sensitivity of the cost function to the control signal
deltaJMin = 1; % min change in cost
currTime = 0; % t_curr
predictionHorizon = endTime / 10; % T

totalN = (endTime / predictionHorizon) * N;
dtval = endTime / totalN;

defaultControlDuration = predictionHorizon/5; % deltaT_init
numItersControlDurationDefault = round(defaultControlDuration*N/predictionHorizon);
nominalControlU1 = 0; % nominal control u1
omega = 0.5; % scale factor omega
samplingTime = 1; % t_s
numIndstCalcMax = 2;
tCalcMax = dtval*numIndstCalcMax; % t_calc max time for iterative control calculations
maxBacktrackingIters = 1; % k_max

Q = [100 0 0; 0 15 0; 0 0 1];
R = [0.1 0; 0 0.1];
P1 = [5 0 0;0 1 0; 0 0 0.01];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initCondsVect = [0 0 pi/2 1 -0.5];
u1init = 1;
u2init = -0.5;

for i = 1:totalN
    u1col(i) = u1init;
    u2col(i) = u2init;
    timevals(i) = dtval * i;
    xdest(:,i) = Fdest(i);
end %% setting init controls

x0col(1) = initCondsVect(1);
x1col(1) = initCondsVect(2);
theta0col(1) = initCondsVect(3);
%%%% setting x init
for i = 2:totalN
    xvalToUseForInitVectDot = [x0col(i-1), x1col(i-1), theta0col(i-1)];
    xvalsVector = Fsinglevectdot(xvalToUseForInitVectDot, [u1init u2init]);
    x0col(i) = x0col(i-1) + dtval * xvalsVector(1);
    x1col(i) = x1col(i-1) + dtval * xvalsVector(2);
    theta0col(i) = theta0col(i-1) + dtval * xvalsVector(3);
end %% setting init trajectory

allControlsInit = [u1col; u2col];
global allPosInit;
allPosInit = [x0col; x1col; theta0col];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% want to find a control vector in each time window
%%%  then find the time tau at which to apply that control vector
%%%  then find the application duration lambda

xsFin = [initCondsVect(1); initCondsVect(2); initCondsVect(3)];
global rhos0;
rhos0 = P1;
rhos(1,:) = [P1(1,:) P1(2,:) P1(3,:)];
xs = [];
window = 1;

while currTime < endTime %%%% requires that endtime is int multiple of pred horizon
    t0 = currTime
    tf = currTime + predictionHorizon;
    xinit = xsFin(:,end);
    % now get a matrix of xs for times t0->tf
    xwindow = getXwindow(xinit, t0, tf, N);
    xs = [xs; xwindow];
    xsFin = transpose(xs(end,:));
    len = length(xs)
    %
    % now get rho (same as adjoint variable P) values for times t0->tf
    rhowindow = getPwindow(xwindow, t0, tf, N);
    rhos = [rhowindow; rhos];
    rhos0 = [rhos(1,1:3); rhos(1,4:6); rhos(1,7:9)];    
    %
    % now get an initial cost J1, init
    tmp = [transpose(xwindow); allControlsInit(:,1:N)];
    J1init = J(tmp)
    %
    % specify a sensitivity alphaD
    %alphaD = alphaD;
    %
    % calculate gamma
    hx = h(xwindow);
    gammas = getGamma(hx, rhowindow);
    %
    % get u_2_*
    u2star = getu2star(gammas, xwindow, hx, rhowindow);
    %
    % now specify tau application time tau > t0 + tcalc
    additionalTauInd = 2;
    tau = t0 + tCalcMax + additionalTauInd*dtval; % some tau value
    %
    % now saturate u2star
    u2star = saturate(u2star, numIndstCalcMax + additionalTauInd);
    %
    currTime = tf;
end



rhos = transpose(rhos);
linspace(0,endTime, totalN);
rhos(1,:);
plot([0:length(rhos(1,:))-1], [rhos(1,:); rhos(2,:); rhos(3,:); rhos(4,:); rhos(5,:); rhos(6,:); rhos(7,:); rhos(8,:); rhos(9,:)]);
xlim([0 length(rhos(1,:))-1]); 
ylim([-20 20]);
title("Position");
xlabel("t");% [allPosInit(1,:);allPosInit(2,:);allPosInit(3,:)]
ylabel("rho vals");
% xs(:,1)
% plot([0:length(xs(:,1))-1], [transpose(xs(:,1)); transpose(xs(:,2)); transpose(xs(:,3))]);
% xlim([0 length(xs(:,1))-1]); 
% ylim([-20 20]);
% title("Position");
% xlabel("t");% [allPosInit(1,:);allPosInit(2,:);allPosInit(3,:)]
% ylabel("rho vals");



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

function hval = h(xs)
  hs = [];
  for i=1:length(xs)
      thetaval = xs(i,3);
      hs(:,:,i) = [cos(thetaval) 0; sin(thetaval) 0; 0 1];
  end
  hval = hs;
end

function Amatv = Amat(x, u, index)
  % D_1(f(x,u))
  Amatv = [0 0 -sin(x(3))*u(1); 0 0 cos(x(3))*u(1); 0 0 0];
end

function xvectdot = Fsinglevectdot(x, u)
  xvectdot = [cos(x(3)) * u(1); sin(x(3)) * u(1); u(2)];
end

function sat = saturate(us, tau)
  global maxU1 maxU2 minU1 minU2 numItersControlDurationDefault
  newUs = us
  for i=tau:(tau+numItersControlDurationDefault)
      if us(i,1) > maxU1
          newU(i,1) = maxU1;
      end
      if us(i,2) > maxU2
          newU(i,2) = maxU2;
      end
      if us(i,1) < minU1
          newU(i,1) = minU1;
      end
      if us(i,2)< minU2
          newU(i,2) = minU2;
      end
  end
  sat = newU;
end

function xs = getXwindow(xinit, t0, tf, N)
  dtval = (tf-t0)/N;
  times = linspace(t0, tf, N);
  xs(1,:) = xinit;
  uinit = [1; -0.5];
  for i=2:N
      xi = xs(i-1,:);
      xs(i,:) = [xi(1) + dtval * cos(xi(3)) * uinit(1); xi(2) + dtval * sin(xi(3)) * uinit(1); xi(3) + dtval * uinit(2)];
  end
end

function gs = getGamma(hx, rhowindow)
  global N;
  gs = [];
  for i=1:N
      hxi = hx(:,:,i);
      rhoi = rhowindow(i,:);
      rhoToUse = [rhoi(1,1:3);rhoi(1,4:6);rhoi(1,7:9)];
      gs(:,:,i) = transpose(hxi)*rhoToUse*transpose(rhoToUse)*hxi;
  end
end

function u2star = getu2star(gammas, xwindow, hs, rhowindow)
  global alphaD N R;
  uinit = [1; -0.5];
  us = [];
  RT = transpose(R);
  xts = transpose(xwindow);
  for i=1:N alphaD;
      gammai = gammas(:,:,i);
      xti = xts(:,i);
      hi = hs(:,:,i);
      rhoi = rhowindow(i,:);
      rhoToUse = [rhoi(1,1:3);rhoi(1,4:6);rhoi(1,7:9)];
      %%%%%%%%%%%%% had to make a vector out of alphaD to make math work
      alphaDVect = [alphaD; alphaD; alphaD];
      %%%%%%%%%%%%%
      us(i,:) = inv(gammai + RT) * (gammai * uinit + transpose(hi) * rhoToUse * alphaDVect);
  end
  u2star = us;
end
%%%
function ps = getPwindow(xwindow, t0, tf, N)
  global allPosInit P1 rhos0
  dtval = (tf-t0)/N;
  times = linspace(t0, tf, N);
  uinit = [1; -0.5];
  rhos = [];
  As = [];
  xidiffs = [];
  for i=1:N
      xi = xwindow(i,:);
      xid = allPosInit(:,i);
      xid;
      xi;
      xidiffs(i,:) = transpose(xi) - xid;
      As(:,:,i) = Amat(xi, uinit, i);
  end
  [Trho rhos] = ode45(@(t,rhos)getRho(t,rhos,As,xidiffs, t0, tf), linspace(tf, t0, N), rhos0);
  ps = flip(rhos);
end
%%%
function rhodot = getRho(t, rho, As, xdiffs, t0, tf)
  global Q N;
  index = round(((t-t0)/(tf-t0))*(N-1)+1);
  rho(1:3);
  rhoToUse = [transpose(rho(1:3));transpose(rho(4:6));transpose(rho(7:9))];
  rhodot = -transpose(As(:,:,index))*rhoToUse - Q * transpose(xdiffs(index,:));
  rhodot = rhodot(:);
end