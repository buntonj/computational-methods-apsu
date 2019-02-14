A = importdata('ans.txt');
info = importdata('info.txt');
n = info(1);
tsteps = info(2);
dx = info(3);
div = info(4);
%T = zeros(n,n,n,tsteps/div);
T = zeros(n,n,n,1);

count = 1;
%for time = 1:tsteps/div
    for i = 1:n
        for j = 1:n
            for k = 1:n
                T(A(count,1),A(count,2),A(count,3),A(count,4)) = A(count,5);
                count = count+1;
            end
        end
    end
%end

limits = [0;0.1;0;0.1;30;160];
[x y] = meshgrid(0:dx:0.1+dx);
mid = ceil(n/2);
%for i = 1:tsteps/div
%surf(x,y,squeeze(T(:,:,mid,i)));
surf(x,y,squeeze(T(:,n,:,1)));
%axis(limits)
drawnow
%end
xlabel('j'); 
ylabel('i'); 
zlabel('Temperature (K)');
title('Temperature at Centered k-coordinate');