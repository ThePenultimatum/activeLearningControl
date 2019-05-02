clear
% initializations 
global epsilon;
epsilon = 0.1;

%eta0 = [[0 0 pi/2 1 -0.5]];

initCondsVect = [0 0 pi/2 1 -0.5];

%j0 = cost of initial eta0

alpha = [0 0.5];
beta = [0 1];

%i = 0;

global Q;
global R;
global P1;
global N;
global T;
global dtval;

Q = [1 0 0; 0 1 0; 0 0 1];
R = [0.1 0; 0 0.1];
P1 = [1 0 0;0 1 0; 0 0 0.01];

N = 100;
T = 2*pi;
dtval = T / N;

%eta0 = [[0 0 0 0 0]];

gamma0 = [[1 1 1 1 1]];

gamma = gamma0;
normGamma = normOfRowVect(gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u1init = 1;
u2init = -0.5;

timevals = [];
%%%%%%%%%%%%%%%%%%%%% setting u init
for i = 1:N
    u1col(i) = u1init;
    u2col(i) = u2init;
    timevals(i) = dtval * i;
end

x0col(1) = initCondsVect(1);
x1col(1) = initCondsVect(2);
theta0col(1) = initCondsVect(3);
%%%% setting x init
for i = 2:N
    xvalToUseForInitVectDot = [x0col(i-1), x1col(i-1), theta0col(i-1)];
    xvalsVector = Fsinglevectdot(xvalToUseForInitVectDot, [u1init u2init]);
    x0col(i) = x0col(i-1) + dtval * xvalsVector(1);
    x1col(i) = x1col(i-1) + dtval * xvalsVector(2);
    theta0col(i) = theta0col(i-1) + dtval * xvalsVector(3);
end

allControlsInit = [u1col; u2col];
allPosInit = [x0col; x1col; theta0col];

initialTrajectory = [allPosInit; allControlsInit];

xcolsToUse = allPosInit;
ucolsToUse = allControlsInit;
xdotcolsToUse = [];

xdestcolsToUse = allPosInit;

Amats = [];
Bmats = [];
amats = [];
bmats = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





while normGamma > epsilon
    for i=1:N
      xtvalsToUse = xcolsToUse(:,i)
      utvalsToUse = ucolsToUse(:,i)
      xtdotcolsToUse = Fvectdot(xcolsToUse, ucolsToUse, i);
      xtdotcolsToUse;
      xdotcolsToUse(:,i) = xtdotcolsToUse;
      At = Amat(xcolsToUse, ucolsToUse, i);
      Amats(:,:,i) = At;
      Bt = Bmat(xcolsToUse, ucolsToUse, i);
      Bmats(:,:,i) = Bt;
%       aval = (transpose(xtvalsToUse-xdestcolsToUse) * Q);
%       amats(:,:,i) = aval;
%       bval = (transpose(utvalsToUse) * R);
%       bmats(:,:,i) = bval;
    end
    [TP, P] = ode45(@(t,P)solvepval(t, P, Q, R, xcolsToUse, ucolsToUse), [10, 0], P1)
    %ode45(@(t,Ps)solvepvalMats(t, Ps, Amats, Bmats, Q, R, xcolsToUse, ucolsToUse), [10, 0], P1)
    %end
    %for i=1:N
    %  [TP, Pt] = ode45(@(t,Pt)solvepval(t, Pt, At, Bt, Q, R), [10, 0], P1)
    %  [Tr, rt] = ode45(@(t,rt)solverval(t, P, At, Bt, Q, R, rt, aval, bval), [10, 0], P1)
    %end  
    normGamma = normGamma - 1;
end

amats
bmats









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function j = J(vals)
  % vals in form of:
  % [[x10, x11, ...., x1N-1],
  %  [x20, x21, ...., x2N-1],
  %  [x30, x31, ...., x3N-1],
  %  [u10, u11, ...., u1N-1],
  %  [u20, u21, ...., u2N-1]]
  xs = vals(1:3, 1:N);
  us = vals(4:5, 1:N);

  global N
  global Q
  global R
  global P1
  Q = [1 0 0; 0 1 0; 0 0 1];
  R = [1 0; 0 1];
  P1 = [10 0 0; 0 10 0; 0 0 10];
  N = 100;
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
  sumsofar = sumsofar + (transpose(xndiff) * P1 * xndiff);
  j = sumsofar;
end

%%%

function dj = DJ(xs, controls, zs, vs)

  global N
  global Q
  global R
  global P1
  Q = [1 0 0; 0 1 0; 0 0 1];
  R = [1 0; 0 1];
  P1 = [10 0 0; 0 10 0; 0 0 10];
  N = 100;
  sumsofar = 0;
  for i = 1:N-1
    xi = [xs(1, i); xs(2, i); xs(3, i)];
    ci = [controls(1, i); controls(2, i)];
    zsi = zs(i,:);
    vsi = vs(i,:);
    sumsofar = sumsofar + transpose(xi) * Q * zsi + transpose(ci) * R * vsi;
  end
  sumsofar;
  sumsofar = sumsofar + transpose(xs(N)) * P1 * zs(N) + 0.5*matrixNorm()
  j = sumsofar
end

%%%

function norm = normOfRowVect(vectVal)
  % this is of the form [[x1,.....,xn]] * transpose([[x1,....,xn]]);
  norm = vectVal * transpose(vectVal);
end

function norm = matrixNorm(matrix, nrows)
  % norm of natrix with n rows
  sumval = 0;
  for i = 0:nrows
      sumval = sumval + normOfRowVect(matrix(i));
  end
  norm = sumval;
end

%%%

function pdot = solvepval(t, P, Q, R, xs, us)
  P
  xs(1,:)
  interp1([1:100], xs(1,:), t)
  interp1([1:100], xs(2,:), t)
  interp1([1:100], xs(3,:), t)
  newx3 = interp1([1:100], xs(3,:), t)
  newu1 = interp1([1:100], us(1,:), t)
  newu2 = interp1([1:100], us(2,:), t)
  A = [0 0 -sin(newx3)*newu1; 0 0 cos(newx3)*newu1; 0 0 0];
  B = [cos(newx3) 0; sin(newx3) 0; 0 1];
  P = reshape(P, size(A));
  pdot = -1*transpose(A)*P - P*A + P*B*(inv(R))*(transpose(B))*P - Q;
  pdot = pdot(:);
end

function pdot = solvepvalMats(t, P, As, Bs, Q, R, xs, us)
  %%%global N;
  %%%N = 100;
  %%%pdot = [];
  
  %for i=1:N
  %  useA = As(:,:,i)
  %  useB = Bs(:,:,i)
  %  P = reshape(P, size(useA));
  %  P;
  %  transpose(useA)*P;
  %  P*useA;
  %  P*useB;
  %  P*useB*(inv(R));
  %  P*useB*(inv(R))*(transpose(useB));
  %  P*useB*(inv(R))*(transpose(useB))*P;
  %  Q;
  %  pdot(:,:,i) = -1*transpose(useA)*P - P*useA + P*useB*(inv(R))*(transpose(useB))*P - Q;
  %end
  
  
  pdot
%pdot = -1*transpose(A)*P - P*A + P*B*(inv(R))*(transpose(B))*P - Q;
pdot = pdot(:);
end

function rdot = solverval(t, P, A, B, R, r, a, b, xs, us)
  rdot = -1*transpose(A)*r - a + P * B * (inv(R)) * b;
  rdot = rdot(:);
end

%%%

function Amatv = Amat(x, u, index)
  % D_1(f(x,u))
  Amatv = [0 0 -sin(x(3,index))*u(1,index); 0 0 cos(x(3,index))*u(1,index); 0 0 0];
end

function Bmatv = Bmat(x, u, index)
  % D_2(f(x,u))
  Bmatv = [cos(x(3,index)) 0; sin(x(3,index)) 0; 0 1];
end

%%%

function xvectdot = Fvectdot(x, u, index)
  xvectdot = [cos(x(3, index)) * u(1, index); sin(x(3, index)) * u(1, index); u(2, index)];
end

function xvectdot = Fsinglevectdot(x, u)
  xvectdot = [cos(x(3)) * u(1); sin(x(3)) * u(1); u(2)];
end

function xdest = Fdest(i)
  global dtval
  global N
  global T
  N = 100;
  T = 2 * pi;
  dtval = T / N;

  xdest = [dtval * 2 * i / pi; 0; pi/2];
end