sinterms = importdata('sinvectors.txt');
costerms = importdata('cosvectors.txt');
%eigen = importdata('eigendata.txt');
%info = importdata('info.txt');
%before = importdata('before.txt');
%potential = importdata('potential.txt');
%nmax = info(1);
%depth = info(2);
%sd = info(3);
%l = info(4);
resolution = info(5);
figure()

%plot(costerms(1:resolution,4),costerms(1:resolution,5))
hold on
%plot(costerms(resolution+1:2*resolution,4),costerms(resolution+1:2*resolution,5))
plot(sinterms(1:resolution,4),sinterms(1:resolution,5));
plot(sinterms(101:200,4),sinterms(101:200,5));
plot(sinterms(201:300,4),sinterms(201:300,5));
%plot(sinterms(301:400,4),sinterms(301:400,5));
%plot(sinterms(401:500,4),sinterms(401:500,5));
title('Sine Terms'); xlabel('x (meters)');
figure()
hold on
plot(costerms(1:resolution,4),costerms(1:resolution,5));
plot(costerms(101:200,4),costerms(101:200,5));
plot(costerms(201:300,4),costerms(201:300,5));
%plot(costerms(301:400,4),costerms(301:400,5));
%plot(costerms(401:500,4),costerms(401:500,5));
title('Cosine Terms'); xlabel('x (meters)');
%plot(potential(:,1),potential(:,2))

