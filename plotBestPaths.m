close all;clear;

datafile1 = 'data1.csv'; datafile2 = 'data2.csv';
flag1 = 167; flag2 = 306;

BestDistancePaths1 = 'BestDistancePaths1.csv';
BestDistancePaths2 = 'BestDistancePaths2.csv';
BestDNPaths1 = 'BestDNPaths1.csv';
BestDNPaths2 = 'BestDNPaths2.csv';
BestPaths1 = 'BestPaths1.csv';
BestPaths2 = 'BestPaths2.csv';

d11 = plot_result(datafile1, BestDistancePaths1, flag1,'Answer1 for data1',1);
d12 = plot_result(datafile2, BestDistancePaths2, flag2,'Answer1 for data2',2);

d21 = plot_result(datafile1, BestDNPaths1, flag1,'Answer2 for data1',3);
d22 = plot_result(datafile2, BestDNPaths2, flag2,'Answer2 for data2',4);

d31 = plot_result(datafile1, BestPaths1, flag1,'Answer3 for data1',5);
d32 = plot_result(datafile2, BestPaths2, flag2,'Answer3 for data2',6);

function d = plot_result(datafile, pathfile, flag, mytitle,i)

data = csvread(datafile); 
path = csvread(pathfile);

A = data(1,:);
B = data(end,:);
dataP = sortrows(data(2:end-1,:),5);
figure;
scatter3(A(2),A(3),A(4),'r','o','filled'); %A
hold on;
scatter3(B(2),B(3),B(4),'r','o','filled'); %B
hold on;
scatter3(dataP(1:flag,2),dataP(1:flag,3),dataP(1:flag,4),'.','m');
hold on;
scatter3(dataP(flag+1:end,2),dataP(flag+1:end,3),dataP(flag+1:end,4),'.','b');
hold on;
plot3([A(2),B(2)],[A(3),B(3)],[A(4),B(4)],'k--');
hold on;
d = [];
 for i = 3: path(1,end)+1
    hold on;
    % calculate distance error
    d(i-2) = 0.001*sqrt( (data(path(1,i+1),2) - data(path(1,i),2))^2 ...
            +(data(path(1,i+1),3) - data(path(1,i),3))^2 ...
            +(data(path(1,i+1),4) - data(path(1,i),4))^2 );
    plot3([data(path(1,i),2),data(path(1,i+1),2)],...
       [data(path(1,i),3),data(path(1,i+1),3)],...
       [data(path(1,i),4),data(path(1,i+1),4)],'k'); 

 end

text(A(2),A(3),A(4),'  A', 'fontsize', 14);
text(B(2),B(3),B(4),'  B', 'fontsize', 14);
xlabel('x');
ylabel('y');
zlabel('z');
title(mytitle);

%saveas(gcf,strcat(mytitle,'.eps'))
end



