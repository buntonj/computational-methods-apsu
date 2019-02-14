A = importdata('pos.txt');
B = importdata('mass.txt');
theta =  0:pi/50:2*pi;
xcirc = 6.371*10^6 * cos(theta);
ycirc = 6.371*10^6 * sin(theta);
plot(xcirc,ycirc)
hold on
plot(A(:,1),A(:,2))