
N = 100;
T = 2 * pi;
dtval = T / N;
x0 = 0;
y0 = 0;
theta0 = pi/2;

u1init = 1;
u2init = -0.5;

Q = [1 0 0; 0 1 0; 0 0 1];
R = [1 0; 0 1];
P1 = [10 0; 0 10];
timevals = []
%%%%%%%%%%%%%%%%%%%%%
for i = 1:N
    u1col(i) = u1init;
    u2col(i) = u2init;
    timevals(i) = dtval * i
end

x0col(1) = x0;
x1col(1) = y0;
theta0col(1) = theta0;
for i = 2:N
    xvalToUseForInitVectDot = [x0col(i-1), y0col(i-1), theta0col(i-1)];
    xvalsVector = Fvectdot(xvalToUseForInitVectDot, [u1init u2init]);
    x0col(i) = x0col(i-1) + dtval * xvalsVector(1);
    y0col(i) = y0col(i-1) + dtval * xvalsVector(2);
    theta0col(i) = theta0col(i-1) + dtval * xvalsVector(3);
end
%testval = [transpose(x0col) transpose(y0col) transpose(theta0col)]
allControlsInit = [u1col; u2col];
allPosInit = [x0col; y0col; theta0col];

allPosInit(1,1:100);
allPosInit(2,1:100);
allPosInit(3,1:100);


%figure(1);
%plot(x0col,y0col);
%title("Initial trajectory");
%figure(1);
%plot(timevals, [u1col; u2col]);
%xlim([0 10]);
%ylim([-1 10]);
%title("Initial controls");
%xlabel("time");
%ylabel("control value");

%%%%%%%%%%%%%%%%%%%%%%
%dtval;
%Fdest(1);
%Fvectdot([0; 0; 0], [0; 0]);
%testval = [[allPosInit]; [allControlsInit]]
%testval(1:3, 1:100)
%JCalc1([[allPosInit]; [allControlsInit]])

optvals = fmincon(@JCalc1, [[allPosInit]; [allControlsInit]], [], [], [], [], [], [], @constraints1)

%%%%%%%%%%%%%%%%%%%%%%

%figure(1);
%plot(optvals(1,:),optvals(2,:));
%title("problem 1 trajectory");
%figure(1);
plot(timevals, [optvals(4,:); optvals(5,:)]);
xlim([0 10]);
ylim([-10 10]);
title("opt controls");
xlabel("time");
ylabel("control vals");







%%%%%%%%%%%%%%%%%%%%%%






function xdest = Fdest(i)
global dtval
global N
global T
N = 100;
T = 2 * pi;
dtval = T / N;

xdest = [dtval * 2 * i / pi; 0; pi/2];
end


function xvectdot = Fvectdot(x, u)
xvectdot = [cos(x(3)) * u(1); sin(x(3)) * u(1); u(2)];
end

function j = JCalc1(vals)
xs = vals(1:3, 1:100);
us = vals(4:5, 1:100);

global N
global Q
global R
global P1
Q = [1 0 0; 0 1 0; 0 0 1];
R = [1 0; 0 1];
P1 = [10 0 0; 0 10 0; 0 0 10];
N = 100;
sumsofar = 0
for i = 1:N-1
    xi = [xs(1, i); xs(2, i); xs(3, i)];
    ui = [us(1, i); us(2, i)];
    xdi = Fdest(i);
    xidiff = xi - xdi;
    sumsofar = sumsofar + transpose(xidiff) * Q * xidiff + transpose(ui) * R * ui;
end
sumsofar
xsn = [xs(1, N); xs(2, N); xs(3, N)]
xdn = Fdest(N)
xndiff = xsn - xdn
sumsofar = sumsofar + (transpose(xndiff) * P1 * xndiff);
j = sumsofar
end

function [c,ceq] = constraints1(vals)
xs = vals(1:3, 1:100);
us = vals(4:5, 1:100);

c = []

global N
global Q
global R
global P1
Q = [1 0 0; 0 1 0; 0 0 1];
R = [1 0; 0 1];
P1 = [10 0 0; 0 10 0; 0 0 10];
N = 100;
T = 2 * pi;
dtval = T / N;

dest1 = Fdest(1)
ceq = [xs(1, 1)-dest1(1), xs(2, 1)-dest1(2), xs(3,1)-dest1(3)]
for i = 2:N-1
    xival = xs(1,i);
    yival = xs(2,i);
    tival = xs(3,i);
    u1ival = us(1,i);
    u2ival = us(2,i);
    vval = [xival yival tival];
    uvval = [u1ival u2ival];
    vdotval = Fvectdot(vval, uvval);
    ceq = [ceq; xs(1,i)-xs(1,i-1)-vdotval(1)*dtval, xs(2,i)-xs(2,i-1)-vdotval(2)*dtval, xs(3,i)-xs(3,i-1)-vdotval(3)*dtval];
end
destN = Fdest(N);
ceq = [ceq; xs(1,N)-destN(1), xs(2,N)-destN(2), xs(3,N)-destN(3)];
end