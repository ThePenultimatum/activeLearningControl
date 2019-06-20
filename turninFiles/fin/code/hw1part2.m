
N = 100;
T = 10;
dtval = T / N;
x0 = 10;
y0 = 0;

x0vect = [10; 0];

Q = [2 0; 0 0.01];
R = [0.1];
P1 = [1 0; 0 0.01];
A = [0 1; -1.6 -0.4];
B = [0;1];
timevals = [];
%%%%%%%%%%%%%%%%%%%%%


solution = bvp4c(@derivativeFunct, @bounds, bvpinit(linspace(0, T, N), [10 0 0 0]))

timevals = solution.x
solution.y
posvars = solution.y(1:2,:)
pvars = solution.y(3:4,:)
controls = -1*inv(R)*transpose(B)*pvars

%plot(timevals, [posvars(1,:); posvars(2,:)]);
%xlim([0 15]);
%ylim([-50, 50]);
%title("opt position");
%xlabel("time");
%ylabel("position vals");


Aperturbs = []
Bperturbs = []
Cperturbs = []
Dperturbs = []
newVectors = []

for i=1:10
    Aperturbs(i) = rand
    Bperturbs(i) = rand
    Cperturbs(i) = rand
    Dperturbs(i) = rand
    newv(i,:) = Aperturbs(i) * sin(timevals * Bperturbs(i) + Cperturbs(i)) + Dperturbs(i);
end

newv

% now have to solve for z and then apply the perturbations in the solution
% to get the derivative parts for the table using newv from abov which was
% calculated with the random perturbations

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
[A -1*(B/R)*(B'); -Q -A']
ts
vect
[A -1*(B/R)*(B')]
[A -1*(B/R)*(B'); -Q -A']
[-Q -A']
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


function j = JCalc2(xs, controls, zs, vs)

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
sumsofar = sumsofar + transpose(xs(N)) * P1 * zs(N)
j = sumsofar
end