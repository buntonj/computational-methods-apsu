A = importdata('ans.txt');
dim = size(A);
n = int32(dim(1)^(1/3));
dx = 0.1/(dim(1)^(1/3)-1);


P = zeros(n,n,n);
for i = 1:n
for j = 1:n
for k = 1:n
P(i,j,k) = A((i-1)*n*n + (j-1)*n + k, 4);
end
end
end

[x y] = meshgrid(0:dx:0.1);
mid = ceil(n/2);
surf(x,y,P(:,:,mid));
xlabel('j'); 
ylabel('i'); 
zlabel('Temperature (K)');
title('Temperature at Centered k-coordinate');

figure();
surf(x,y,squeeze(P(mid,:,:)));
xlabel('k');
ylabel('j');
zlabel('Temperature (K)');
title('Temperature at Centered i-coordinate');

figure();
surf(x,y,squeeze(P(:,mid,:)));
xlabel('k');
ylabel('i');
zlabel('Temperature (K)');
title('Temperature at Centered j-coordinate');


filename='Tempdistr'
limits = [0;0.1;0;0.1;2;350];
h = figure;
for i = 1:n
%figure()
surf(x,y,squeeze(P(i,:,:)));
title(num2str(i))
axis(limits)
drawnow
% Capture the plot as an image
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
% Write to the GIF File
if i == 1
imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.05);
else
imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.05);
end
end

