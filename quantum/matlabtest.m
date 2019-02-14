elevels = 3;
info = importdata('info.txt');
potential = importdata('potential.txt');
nmax = info(1);
depth = info(2);
sd = info(3);
l = info(4);
resolution = info(5);
hbar = 1.0545718E-34;
m = 9.1094E-31;
x = -l:2*l/resolution:l;

Hsin = zeros(nmax,nmax);
Hcos = zeros(nmax+1,nmax+1);
sinpsi = zeros(resolution,elevels);
cospsi = zeros(resolution,elevels);
fxsin = @(x) depth*sin(n*pi*x/l).*sin(m*pi*x/l).*exp(-sd*(x/l).^2);
fxcos = @(x) depth*cos(n*pi*x/l).*cos(m*pi*x/l).*exp(-sd*(x/l).^2);

for a = 1:nmax
    n = a;
    for b = 1:nmax
        m = b;
        Hsin(a,b) = (2/l)*integral(fxsin,-l,l);
        Hsin(b,a) = Hsin(a,b);
    end
    Hsin(a,a) = (pi*hbar*a/l)^2/(2*m) + Hsin(a,a);
end

[vec val] = eig(Hsin);

for a = 1:elevels
    for b = 1:size(vec,2)
        for c = 1:resolution
        sinpsi(c,a) = sinpsi(c,a) + sqrt(2/l)*vec(b,a)*sin(b*pi*x(c)/l);
        end
    end
    sinpsi(:,a) = sinpsi(:,a)/norm(sinpsi(:,a));
end

Hcos = zeros(nmax+1,nmax+1);

for a = 0:nmax
    n = a;
    for b = 0:nmax
        m = b;
        Hcos(a+1,b+1) = (2/l)*integral(fxcos,-l,l);
        Hcos(b+1,a+1) = Hcos(a+1,b+1);
    end
    Hcos(a+1,a+1) = (pi*hbar*a/l)^2/(2*m) + Hcos(a+1,a+1);
end

[vec val] = eig(Hcos);

for a = 1:elevels
    for b = 1:size(vec,2)
        for c = 1:resolution
        cospsi(c,a) = cospsi(c,a) + sqrt(2/l)*vec(b,a)*cos((b-1)*pi*x(c)/l);
        end
    end
    cospsi(:,a) = cospsi(:,a)/norm(cospsi(:,a));
end

figure();
plot(x(1:size(x,2)-1),sinpsi(:,1))
%hold on
%plot(x(1:size(x,2)-1),cospsi(:,1))
%plot(potential(:,1),potential(:,2))
title('Matlab Result')   