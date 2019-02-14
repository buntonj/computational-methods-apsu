A = importdata('sinsinsinvectors.txt');
info = importdata('info.txt');
nmax = info(1);
depth = info(2);
sd = info(3);
l = info(4);
resolution = info(5);
elevels = 3;

sinsinsin = zeros(resolution+1,resolution+1,resolution+1,elevels);

count = 1;
for elevel = 1:elevels
for i = 1:resolution+1
    for j = 1:resolution+1
        for k = 1:resolution+1
            sinsinsin(A(count,2)+1,A(count,3)+1,A(count,4)+1,elevels) = A(count,5);
        end
    end
end
end

            