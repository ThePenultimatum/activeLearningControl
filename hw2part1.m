clear
% initializations 
global Q R P1 N T dtval epsilon invR allPosInit;

epsilon = 0.1;

initCondsVect = [0 0 pi/2 1 -0.5];

alpha = 0.4;
beta = 0.7;

Q = [100 0 0; 0 15 0; 0 0 1];
R = [0.1 0; 0 0.1];
invR = inv(R);
P1 = [5 0 0;0 1 0; 0 0 0.01];

N = 1000;
T = 2*pi;
dtval = T / N;

norm = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u1init = 1;
u2init = -0.5;

global xdest;

%%%%%%%%%%%%%%%%%%%%% setting u init
u1col = zeros(1,N);
u2col = zeros(1,N);
timevals = zeros(1,N);
xdest = zeros(3,N);
for i = 1:N
    u1col(i) = u1init;
    u2col(i) = u2init;
    timevals(i) = dtval * i;
    xdest(:,i) = Fdest(i);
end
x0col = zeros(1,N);
x1col = zeros(1,N);
theta0col = zeros(1,N);
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

Amats = zeros(3,3,N);
Bmats = zeros(3,2,N);
amats = [];
bmats = [];

z0 = [0; 0; 0];
v0 = [0; 0];
zs = [];
vs = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iters = 0;
while (abs(norm) > epsilon) && (iters < 100)
    
    %%% calculate 
    for i=1:N
      At = Amat(xcolsToUse, ucolsToUse, i);
      Amats(:,:,i) = At;
      Bt = Bmat(xcolsToUse, i);
      Bmats(:,:,i) = Bt;
    end
    
    [TP, P] = ode45(@(t,P)solvepval(t, P, Q, Amats, Bmats), linspace(T,0,N), P1);
    tmpP = P(N,:);
    
    r0 = [tmpP(1,1:3); tmpP(1,4:6); tmpP(1, 7:9)] * xcolsToUse(:,length(xcolsToUse));%z0; %P1 * (xcolsToUse(:,N) - xdest(:,N));%[tmpP(1,1:3); tmpP(1,4:6); tmpP(1, 7:9)] * z0;
    [Tr, r] = ode45(@(t,r)solverval(t, r, P, R, Q, Amats, Bmats, xcolsToUse, ucolsToUse), linspace(T,0,N), r0);
    
    x0valForZdot = [0; 0; 0];%[xcolsToUse(1,1); xcolsToUse(2,1); xcolsToUse(3,1)];
    [Tz, z] = ode45(@(t,z)getZdot(t, z, flip(P), R, Amats, Bmats, ucolsToUse, flip(r)), linspace(0,T,N), x0valForZdot);

    v = getV(R, Bmats, flip(P), z, flip(r), ucolsToUse);
    

    gamma = 0.001;
    newucols = ucolsToUse + v * gamma;
    newxcols = zeros(3,N);
    prev = allPosInit(:,1);
    newxcols(:,1) = [prev(1); prev(2); prev(3)];
    for i=2:(T/dtval)
        dotControls = newucols(:,i-1);
        ustousehere = Fsinglevectdot(prev, dotControls);
        newx1 = prev(1) + dtval * ustousehere(1);
        newx2 = prev(2) + dtval * ustousehere(2);
        newx3 = prev(3) + dtval * ustousehere(3);
        newxcols(:,i+1) = [newx1; newx2; newx3];
        prev = [newx1; newx2; newx3];
    end
    
    
    xcolsToUse = newxcols;
    ucolsToUse = newucols;
    
    %%%
    
    norm = DJ(xcolsToUse, ucolsToUse, z, v) %matrixNorm([transpose(z); v], 5)
    iters = iters + 1
end

% amats;
% bmats;
% length(P);
% P = transpose(P);
% xcolsToUse(2,:);
% plot([0:N-1], [P(1,:); P(2,:); P(3,:); P(4,:); P(5,:); P(6,:); P(7,:); P(8,:); P(9,:)]); %[xcolsToUse(1,:);xcolsToUse(2,:);xcolsToUse(3,:)]);%
% xlim([0 N-1]); 
% ylim([-100 100]);
% title("p");
% xlabel("time");% [allPosInit(1,:);allPosInit(2,:);allPosInit(3,:)]
% ylabel("p vals");

% plot(dtval * [0:N-1],[ucolsToUse(1,:);ucolsToUse(2,:)]);
% xlim([0 dtval*(N-1)]); 
% ylim([-35 35]);
% title("Controls");
% xlabel("time");% [allPosInit(1,:);allPosInit(2,:);allPosInit(3,:)]
% ylabel("Control values");

% plot(dtval * [0:N-1],[xcolsToUse(1,:);xcolsToUse(2,:);xcolsToUse(3,:)]);
% xlim([0 dtval*(N-1)]); 
% ylim([-35 35]);
% title("Position");
% xlabel("time");% [allPosInit(1,:);allPosInit(2,:);allPosInit(3,:)]
% ylabel("Position values");

plot(xcolsToUse(1,:),xcolsToUse(2,:));
xlim([-5 5]); 
ylim([-5 5]);
title("Position");
xlabel("x1");% [allPosInit(1,:);allPosInit(2,:);allPosInit(3,:)]
ylabel("x2");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  xsn = [xs(1, N); xs(2, N); xs(3, N)];
  xdn = Fdest(N);
  xndiff = xsn - xdn;
  sumsofar = 0.5*sumsofar + 0.5*(transpose(xndiff) * P1 * xndiff);
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
    vsi = vs(:,i);
    xdi = Fdest(i);
    xidiff = xi - xdi;
    sumsofar = sumsofar + transpose(xidiff) * Q * transpose(zsi) + transpose(ci) * R * vsi;
  end
  xsn = [xs(1, N); xs(2, N); xs(3, N)];
  xdn = Fdest(N);
  xndiff = xsn - xdn;
  sumsofar = 0.5*sumsofar + 0.5*(transpose(xndiff) * P1 * transpose(zs(N,:)));
  dj = sumsofar;
end

%%%

function pdot = solvepval(t, P, Q, As, Bs)
  global N T invR;
  index = round((t/T)*(N-1)+1);
  A = As(:,:,index);
  B = Bs(:,:,index);
  
  P = reshape(P, size(A));
  pdot = -transpose(A)*P - P*A + P*B*(invR)*(transpose(B))*P - Q;
  pdot = pdot(:);
end

function rdot = solverval(t, r, P, R, Q, As, Bs, xs, us)
  global xdest N T invR;
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
  
  rdot = -1*transpose(A - B * invR * transpose(B) * newP)*r - a + newP * B * (invR) * b;
  rdot = rdot(:);
  
end

function zdot = getZdot(t, z, P, R, As, Bs, us, rs)
  global N T invR;
  index = round((t/T)*(N-1)+1);
  A = As(:,:,index);
  B = Bs(:,:,index);
  newu1 = us(1,index);
  newu2 = us(2,index);
  newr1 = rs(index,1);
  newr2 = rs(index,2);
  newr3 = rs(index,3); 
  b = transpose(transpose([newu1; newu2]) * R);
  Pval = P(index,:);
  newP = [Pval(1:3); Pval(4:6); Pval(7:9)];
  rval = [newr1; newr2; newr3];
  
  zdot = A * z + B * (-1 * invR * transpose(B) * newP * z - invR * transpose(B) * rval - invR * b);
end

function v = getV(R, Bs, P, zs, rs, us)
global invR;
  vs = zeros(2,length(rs));
  for i=1:length(rs)
    B = Bs(:,:,i);
    Pval = P(i,:);
    newP = [Pval(1:3); Pval(4:6); Pval(7:9)];
    newu1 = us(1,i);
    newu2 = us(2,i);
    b = transpose(transpose([newu1; newu2]) * R);
    z = transpose(zs(i,:));
    rval = transpose(rs(i,:));
    vs(:,i) = -invR * transpose(B) * newP * z - invR * transpose(B) * rval - invR * b;
  end
  v = vs;%-1 * invR * transpose(B) * newP * z - invR * transpose(B) * rval - invR * b;
end

%%%

function Amatv = Amat(x, u, index)
  % D_1(f(x,u))
  Amatv = [0 0 -sin(x(3,index))*u(1,index); 0 0 cos(x(3,index))*u(1,index); 0 0 0];
end

function Bmatv = Bmat(x,index)
  % D_2(f(x,u))
  Bmatv = [cos(x(3,index)) 0; sin(x(3,index)) 0; 0 1];
end

function xvectdot = Fsinglevectdot(x, u)
  xvectdot = [cos(x(3)) * u(1); sin(x(3)) * u(1); u(2)];
end

function xdest = Fdest(i)
  global dtval;

  xdest = [dtval * 2 * i / pi; 0; pi/2];
end