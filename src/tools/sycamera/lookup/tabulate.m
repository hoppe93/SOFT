%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE LOOKUP TABLE FOR \Zeta_{mn}
% 
% Written by: Mathias Hoppe (2018)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear;

% PARAMETERS
VISUALIZE=false;
OUTPUT=false;

SAVEPATH='./';

% Physical limits
gammaMax = 200;
kmin = 1e-6;

% Resolution
nx = 100;
ny = 100;
xbreak = 0.9;
ybreakl = 0.04;
ybreaku = 0.98;
nsmoothx = 50;
nsmoothy = 50;
epsy = 1-1e-6;  % How close to 1 we should tabulate a logarithmically distributed value in y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERNAL. DON'T TOUCH.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MaxX = sqrt(gammaMax^2-1) / gammaMax;
MinY = (2*kmin*(gammaMax-sqrt(gammaMax^2-1))) / (((gammaMax-kmin)^2-1) + 2*kmin*(gammaMax-sqrt(gammaMax^2-1)));

% Create grid
if nsmoothx > nx || nsmoothy > ny
    error('The smooth part of the grid cannot have more than the total number of points.');
elseif mod(ny,2) ~= 0 || mod(nsmoothy,2) ~= 0
    error('The y-resolution must be an even number.');
end

lx = linspace(0, xbreak, nsmoothx+1);
ly = linspace(ybreakl, ybreaku, nsmoothy+2);
x = [lx(1:end-1), 1-logspace(log10(1-xbreak), log10(1-MaxX), nx-nsmoothx)];
y = [logspace(log10(MinY), log10(ybreakl), (ny-nsmoothy)/2), ly(2:end-1), logspace(log10(ybreaku), log10(epsy), (ny-nsmoothy)/2-1), 1];

[X,Y] = meshgrid(x, y);

% Evaluate function
I11 = zeros(nx,ny);
I12 = zeros(nx,ny);
I21 = zeros(nx,ny);
for i=1:nx
    for j=1:ny
        I11(i,j) = ZetaMN(1, 1, X(i,j), Y(i,j));
        I12(i,j) = ZetaMN(1, 2, X(i,j), Y(i,j));
        I21(i,j) = ZetaMN(2, 1, X(i,j), Y(i,j));
    end
end

if VISUALIZE
    figure(1);
    surf(X,Y,I11);
    set(gca, 'TickLabelInterpreter', 'latex');
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$y$', 'Interpreter', 'latex');
    title('$\mathrm{Z}_{11}(x,y)$', 'Interpreter', 'latex');
    
    figure(2);
    surf(X,Y,I12);
    set(gca, 'TickLabelInterpreter', 'latex');
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$y$', 'Interpreter', 'latex');
    title('$\mathrm{Z}_{12}(x,y)$', 'Interpreter', 'latex');
    
    figure(3);
    surf(X,Y,I21);
    set(gca, 'TickLabelInterpreter', 'latex');
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$y$', 'Interpreter', 'latex');
    title('$\mathrm{Z}_{21}(x,y)$', 'Interpreter', 'latex');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if OUTPUT
    output2DLookup(SAVEPATH,'Zeta11', x, y, I11, kmin, gammaMax);
    output2DLookup(SAVEPATH,'Zeta12', x, y, I12, kmin, gammaMax);
    output2DLookup(SAVEPATH,'Zeta21', x, y, I21, kmin, gammaMax);
end

disp('Done.');

%% Test interpolation
clc;

% Generate shifted grid
xp = zeros(1,nx-1);
yp = zeros(1,ny-1);
for i=1:nx-1
    xp(i) = (x(i)+x(i+1))/2;
end
for i=1:ny-1
    yp(i) = (y(i)+y(i+1))/2;
end

[Xp, Yp] = meshgrid(xp, yp);

Vq = interp2(X, Y, I11, Xp, Yp, 'makima');
Zp = zeros(nx-1, ny-1);
for i=1:nx-1
    for j=1:ny-1
        Zp(i,j) = ZetaMN(1,1,Xp(i,j),Yp(i,j));
    end
end

err = abs((Vq-Zp)./Vq);

figure(4);
surf(Xp, Yp, err);
xlabel('x');
ylabel('y');
