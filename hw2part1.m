clear
% initializations 
global epsilon;
epsilon = 0.1;

%eta0 = [[0 0 pi/2 1 -0.5]];

initCondsVect = [0 0 pi/2 1 -0.5];

%j0 = cost of initial eta0

alpha = 0.4;
beta = 0.7;

%i = 0;

global Q R P1 N T dtval;

Q = [100 0 0; 0 15 0; 0 0 1];
R = [0.1 0; 0 0.1];
P1 = [5 0 0;0 1 0; 0 0 0.01];

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

global xdest;
xdest = [];

%%%%%%%%%%%%%%%%%%%%% setting u init
for i = 1:N
    u1col(i) = u1init;
    u2col(i) = u2init;
    timevals(i) = dtval * i;
    xdest(:,i) = Fdest(i);
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
global allPosInit;
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

z0 = [0; 0; 0];
v0 = [0; 0];
zs = [];
vs = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


while normGamma > epsilon
    for i=1:N
      xtvalsToUse = xcolsToUse(:,i);
      utvalsToUse = ucolsToUse(:,i);
      xtdotcolsToUse = Fvectdot(xcolsToUse, ucolsToUse, i);
      xtdotcolsToUse;
      xdotcolsToUse(:,i) = xtdotcolsToUse;
      At = Amat(xcolsToUse, ucolsToUse, i);
      Amats(:,:,i) = At;
      Bt = Bmat(xcolsToUse, ucolsToUse, i);
      Bmats(:,:,i) = Bt;
    end
      
    %%%%%%%%%%%%%%%% a and b are defined forwards while p and r are defined
    %%%%%%%%%%%%%%%% backwards in time
  
  
    [TP, P] = ode45(@(t,P)solvepval(t, P, Q, R, Amats, Bmats, xcolsToUse, ucolsToUse), linspace(T,0,N), P1)
    tmpP = P(100,:);
    r0 = [tmpP(1,1:3); tmpP(1,4:6); tmpP(1, 7:9)] * z0;
    
    [Tr, r] = ode45(@(t,r)solverval(t, r, P, R, Q, Amats, Bmats, xcolsToUse, ucolsToUse), linspace(T,0,N), r0)
    x0valForZdot = [xcolsToUse(1,1); xcolsToUse(2,1); xcolsToUse(3,1)];
    [Tz, z] = ode45(@(t,z)getZdot(t, z, P, R, Q, Amats, Bmats, xcolsToUse, ucolsToUse, r), linspace(T,0,N), x0valForZdot)
    v = getV(R, Amats, Bmats, P, z, r, xcolsToUse, ucolsToUse, Q)
    
    %%% armijo
    
    n = 0;
    gammaval = beta ^ n;
    
    fvalleft = J(zeta plus);
    fvalright = J(zeta cur) + alpha * beta * DJ(zeta curr in direction of eta calculated above);
    while (fvalleft > fvalright) && n < 15
        n = n + 1;
        gammaval = beta ^ n;
        fvalleft = ;
        fvalright = ;
    end
    
    %%%
    
    
    normGamma = normGamma - 1;
end

amats;
bmats;
length(P)
P = transpose(P);
plot([0:99], [P(1,:); P(2,:); P(3,:); P(4,:); P(5,:); P(6,:); P(7,:); P(8,:); P(9,:)]);
xlim([0 100]);
ylim([-50 50]);
title("p");
xlabel("time");
ylabel("p vals");

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

  global N Q R P1

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

  global N Q R P1

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

function pdot = solvepval(t, P, Q, R, As, Bs, xs, us)
  global N T;
  index = round((t/T)*(N-1)+1);
  P;
  A = As(:,:,index);
  B = Bs(:,:,index);
  
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

function rdot = solverval(t, r, P, R, Q, As, Bs, xs, us)
  global xdest N T;
  index = round((t/T)*(N-1)+1);
  
  newx1 = xs(1,index);
  newx2 = xs(2,index);
  newx3 = xs(3,index);
  newu1 = us(1,index);
  newu2 = us(2,index);
  
  xdestcols1 = xdest(1,index);
  xdestcols2 = xdest(2,index);
  xdestcols3 = xdest(3,index);
  xdestcolsToUse = [xdestcols1; xdestcols2; xdestcols3];
  
  A = As(:,:,index);
  B = Bs(:,:,index);
  
  a = transpose(transpose([newx1; newx2; newx3]-xdestcolsToUse) * Q);
  b = transpose(transpose([newu1; newu2]) * R);
  Pval = P(index,:);
  newP = [Pval(1:3); Pval(4:6); Pval(7:9)];
  
  rdot = -1*transpose(A - B * inv(R) * transpose(B) * newP)*r - a + newP * B * (inv(R)) * b;
  rdot = rdot(:);
  
end

function zdot = getZdot(t, z, P, R, Q, As, Bs, xs, us, rs)
  global xdest N T;
  index = round((t/T)*(N-1)+1);
  A = As(:,:,index);
  B = Bs(:,:,index);
  newx1 = xs(1,index);
  newx2 = xs(2,index);
  newx3 = xs(3,index);
  newu1 = us(1,index);
  newu2 = us(2,index);
  newr1 = rs(index,1);
  newr2 = rs(index,2);
  newr3 = rs(index,3); 
  xdestcols1 = xdest(1,index);
  xdestcols2 = xdest(2,index);
  xdestcols3 = xdest(3,index);
  xdestcolsToUse = [xdestcols1; xdestcols2; xdestcols3];
  a = transpose(transpose([newx1; newx2; newx3]-xdestcolsToUse) * Q);
  b = transpose(transpose([newu1; newu2]) * R);
  Pval = P(index,:);
  newP = [Pval(1:3); Pval(4:6); Pval(7:9)];
  rval = [newr1; newr2; newr3];
  
  zdot = A * z + B * (-1 * inv(R) * transpose(B) * newP * z - inv(R) * transpose(B) * rval - inv(R) * b);
end

function v = getV(R, As, Bs, P, zs, rs, xs, us, Q)
  global T N xdest;
  vs = [];
  for i=1:length(rs)
    A = As(:,:,i);
    B = Bs(:,:,i);
    Pval = P(i,:);
    newP = [Pval(1:3); Pval(4:6); Pval(7:9)];
    newx1 = xs(1,i);
    newx2 = xs(2,i);
    newx3 = xs(3,i);
    newu1 = us(1,i);
    newu2 = us(2,i);
    xdestcols1 = xdest(1,i);
    xdestcols2 = xdest(2,i);
    xdestcols3 = xdest(3,i);
    xdestcolsToUse = [xdestcols1; xdestcols2; xdestcols3];
    a = transpose(transpose([newx1; newx2; newx3]-xdestcolsToUse) * Q);
    b = transpose(transpose([newu1; newu2]) * R);
    z = transpose(zs(i,:));
    rval = transpose(rs(i,:));
    vs(:,i) = -1 * inv(R) * transpose(B) * newP * z - inv(R) * transpose(B) * rval - inv(R) * b;
  end
  v = vs%-1 * inv(R) * transpose(B) * newP * z - inv(R) * transpose(B) * rval - inv(R) * b;
end

function Amatv = Amat(x, u, index)
  % D_1(f(x,u))
  Amatv = [0 0 -sin(x(3,index))*u(1,index); 0 0 cos(x(3,index))*u(1,index); 0 0 0];
end

function Bmatv = Bmat(x, u, index)
  % D_2(f(x,u))
  Bmatv = [cos(x(3,index)) 0; sin(x(3,index)) 0; 0 1];
end


function xvectdot = Fvectdot(x, u, index)
  xvectdot = [cos(x(3, index)) * u(1, index); sin(x(3, index)) * u(1, index); u(2, index)];
end

function xvectdot = Fsinglevectdot(x, u)
  xvectdot = [cos(x(3)) * u(1); sin(x(3)) * u(1); u(2)];
end

function xdest = Fdest(i)
  global dtval;
  global N;
  global T;
  N = 100;
  T = 2 * pi;
  dtval = T / N;

  xdest = [dtval * 2 * i / pi; 0; pi/2];
end