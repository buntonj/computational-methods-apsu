% IMPORT LAST RUN DETAILS
info = importdata('info.txt');
nmax = info(1);
depth = info(2);
sd = info(3);
l = info(4);
resolution = info(5);
%dx = info(6);
dx = 2E-12;
elevels = 5;
[x y] = meshgrid(-l:dx:l);


% IMPORT SIN/SIN/SIN TERMS
A = importdata('sinsinsinvectors.txt');
sinsinsin = zeros(resolution+1,resolution+1,resolution+1,elevels);
sinsinsinenergies = zeros(3,1);
count = 1;
for elevel = 1:elevels
    sinsinsinenergies(elevel) = A(count,1);
for i = 1:resolution+1
    for j = 1:resolution+1
        for k = 1:resolution+1
            sinsinsin(A(count,2)+1,A(count,3)+1,A(count,4)+1,elevel) = A(count,5);
            count = count+1;
        end
    end
end
end

%IMPORT SIN/SIN/COS TERMS
A = importdata('sinsincosvectors.txt');
sinsincos = zeros(resolution+1,resolution+1,resolution+1,elevels);
sinsincosenergies = zeros(3,1);
count = 1;
for elevel = 1:elevels
    sinsincosenergies(elevel) = A(count,1);
for i = 1:resolution+1
    for j = 1:resolution+1
        for k = 1:resolution+1
            sinsincos(A(count,2)+1,A(count,3)+1,A(count,4)+1,elevel) = A(count,5);
            count = count+1;
        end
    end
end
end

%IMPORT SIN/COS/COS TERMS
A = importdata('sincoscosvectors.txt');
sincoscos = zeros(resolution+1,resolution+1,resolution+1,elevels);
sincoscosenergies = zeros(3,1);
count = 1;
for elevel = 1:elevels
    sincoscosenergies(elevel) = A(count,1);
for i = 1:resolution+1
    for j = 1:resolution+1
        for k = 1:resolution+1
            sincoscos(A(count,2)+1,A(count,3)+1,A(count,4)+1,elevel) = A(count,5);
            count = count+1;
        end
    end
end
end

%IMPORT COS/COS/COS TERMS
A = importdata('coscoscosvectors.txt');
coscoscos = zeros(resolution+1,resolution+1,resolution+1,elevels);
coscoscosenergies = zeros(3,1);
count = 1;
for elevel = 1:elevels
    coscoscosenergies(elevel) = A(count,1);
for i = 1:resolution+1
    for j = 1:resolution+1
        for k = 1:resolution+1
            coscoscos(A(count,2)+1,A(count,3)+1,A(count,4)+1,elevel) = A(count,5);
            count = count+1;
        end
    end
end
end         