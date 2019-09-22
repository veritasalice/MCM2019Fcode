close all;
% [data1, datac1] = data_prep('data1.csv',306);
% [data2, datac2] = data_prep('data2.csv',167);
datafile = 'data1.csv';
pathfile = 'BestDistancePaths2.csv';
flag = 167;
data = csvread(datafile); 
path = csvread(pathfile);

N = data(:,1); %?????
X = data(:,2);
Y = data(:,3);
Z = data(:,4);
T = data(:,5); %?????
L = data(:,6);

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
plot3([A(2),B(2)],[A(3),B(3)],[A(4),B(4)],'b');
hold on;
%------------------------------------------------------------
pathn = [1 579 137 171 12 287 613];
for i = 1: 6 %i = 3: 3+pathn(1,end)-1
   hold on;
   j = pathn(i)
   plot3([data(j,2),data(j+1,2)],...
       [data(j,3),data(j+1,3)],[data(j,4),data(j+1,4)]);
 
end
%-------------------------------------------------------------
text(A(2),A(3),A(4),'  A');
text(B(2),B(3),B(4),'  B');
xlabel('x');
ylabel('y');
zlabel('z');

