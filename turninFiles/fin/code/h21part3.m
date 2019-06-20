clear
N = 100;
T = 10;
dtval = T / N;
x0 = 10;
y0 = 0;

x0vect = [10; 0];

Q = [2 0; 0 0.01];
R = 0.1;
P1 = [1; 0; 0; 0.01];
A = [0 1; -1.6 -0.4];
B = [0;1];
timevals = [];
%%%%%%%%%%%%%%%%%%%%%


%solution = bvp4c(@derivativeFunct, @bounds, bvpinit(linspace(0, T, N), [10 0 0 0]))

%timevals = solution.x
%solution.y
%posvars = solution.y(1:2,:)
%pvars = solution.y(3:4,:)
%controls = -1*inv(R)*transpose(B)*pvars


[TP, P] = ode45(@(t,P)riccati(t, P, A, B, Q, R), [10, 0], P1)

[T, x] = ode45(@(T,x)applyriccati(T,x,TP,P,A,B,R), [0,10], x0vect)

%plot(T, [transpose(x(:,1)); transpose(x(:,2))]);
%xlim([0 15]);
%ylim([-50, 50]);
%title("opt position");
%xlabel("time");
%ylabel("position vals");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%problem 2

Q = [2 0; 0 0.01];
R = [0.1];
P1 = [1 0; 0 0.01];
A = [0 1; -1.6 -0.4];
B = [0;1];
timevals = [];



solution = bvp4c(@derivativeFunct, @bounds, bvpinit(linspace(0, 10, 100), [10 0 0 0]))

timevals = solution.x
solution.y
posvars = solution.y(1:2,:)
pvars = solution.y(3:4,:)
controls = -1*inv(R)*transpose(B)*pvars

%%%
controls3 = []
lenval = length(P(:,1))
for i = 1:length(T)
    controls3(i) = -1*inv(R)*transpose(B)*[P(lenval,1) P(lenval,2); P(lenval,3) P(lenval,4)]*transpose(x(i,:))
end

controls;
controls3;
for i = 1:33
    l = length(controls);
    randval = round(rand * 100);
    if randval <= 0
        randval = 1;
    end
    if randval == 100
        randval = 99
    end    
    newval = controls(randval) + controls(randval+1) / 2;
    controls = [controls(1:randval), newval, controls(randval+1:l)];
end
controls;
controls3;
ctrldiff = controls - controls3

posvars = transpose(posvars)
x
length(posvars);
length(x);
for i = 1:33
    l = length(posvars);
    randval = round(rand * 100);
    if randval <= 0
        randval = 1;
    end
    if randval == 100
        randval = 99
    end
    
    newx1 = posvars(randval,1) + posvars(randval+1,1) / 2;
    newx2 = posvars(randval,2) + posvars(randval+1,2) / 2;
    posvars = [posvars(1:randval,:); newx1 newx2 ; posvars(randval+1:l,:)];
end
posdiff = posvars - x

%%%

%plot(T, [transpose(posdiff(:,1)); transpose(posdiff(:,2))]);
%xlim([0 10]);
%ylim([-30, 30]);
%title("opt position differences");
%xlabel("time");
%ylabel("position difference vals");

plot(T, ctrldiff);
xlim([0 10]);
ylim([-30, 30]);
title("opt control differences");
xlabel("time");
ylabel("control difference vals");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









function res = applyriccati(T, x, TP, P, A, B, R)
lenval = length(P(:,1))
res = (A - B * inv(R) * transpose(B) * [P(lenval,1) P(lenval,2); P(lenval,3) P(lenval,4)]) * x
end

function pdot = riccati(t, P, A, B, Q, R)
P;
A;
P = reshape(P, size(A));
pdot = -1*transpose(A)*P - P*A + P*B*(transpose(B))*P - Q;
pdot = pdot(:);
end

function g = guess(x)
g = [sin(x)
     cos(x)];
end
%%%%%%%%%%%%%%%%%%%%%

function xvectdot = Fvectdot(x, u)
global A;
A = [0 1; -1.6 -0.4];
global B;
Q = [0 1];
xvectdot = A * x + B * u;
end

function pvectdot = Fpvectdot(p, x)
global A;
A = [0 1; -1.6 -0.4];
global Q;
Q = [2 0; 0 0.01];
pvectdot = -1 * transpose(A) * p - Q * x;
end

function d = derivativeFunct(ts, vect)

Q = [2 0; 0 0.01];
R = [0.1];
P1 = [1 0; 0 0.01];
A = [0 1; -1.6 -0.4];
B = [0;1];
[A -1*(B/R)*(B'); -Q -A'];
ts;
vect;
d = [A -1*(B/R)*(B'); -Q -A'] * vect;
end

function b = bounds(xs0, xsf)

N = 100;
T = 10;
dtval = T / N;
x0 = 10;
y0 = 0;

Q = [2 0; 0 0.01];
R = [0.1];
P1 = [1 0; 0 0.01];
A = [0 1; -1.6 -0.4];
B = [0;1];
timevals = [];


x0vect = [10; 0];
ps = xsf(3:4)
xvals = xsf(1:2)
ppart = ps - P1 * xvals;
b = [xs0(1)-x0vect(1); xs0(2)-x0vect(2); ppart(1); ppart(2)];
end