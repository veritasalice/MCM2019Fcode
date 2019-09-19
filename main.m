clear;close all;

[data1, datac1] = data_prep('data1.csv',306);
[data2, datac2] = data_prep('data2.csv',167);
 
% % save data
% writematrix(data1, 'data1.csv');
% writematrix(datac1, 'datac1.csv');
% writematrix(data2, 'data2.csv');
% writematrix(datac2, 'datac2.csv');

data = datac2;
n = length(data);
graph = zeros(n,n);

for i = 1:n
    for j = 1:n       
        % i Vertical：1  or Horizontal：0 
        
        % if ix < jx then calculate(j in front of i) 
        if data(i,2) < data(j,2)
            
        
        if data(i,5) == 10 %start A===============================================
            
            % start is A, end is B-----------------------------------------
            if data(j,5) == 100
                graph(i,j) = 1; % save to graph
                
            % start is A, end is not B-------------------------------------
            else
                d = sqrt((data(j,2))^2+(data(j,3))^2+(data(j,4))^2);
                
                if data(j,5) == 1 % vertical                   
                    gamma = 1.5e4;
                    if d > gamma
                        graph(i,j) = 1; % save to graph
                    end

                else % horizontal
                    gamma = 2e4;
                    if d > gamma
                        graph(i,j) = 1; % save to graph
                    end

                end
            end
            
        else %start not A=========================================================
            
            d = sqrt((data(j,2)-data(i,2))^2+(data(j,3)-data(i,3))^2+(data(j,4)-data(i,4))^2);          
            % start is not A, end is B--------------------------------------
            if data(j,5) == 100
                gamma =3e4;
                if  data(j,5) == 1 % vertical
                    if d > gamma
                        graph(i,j) = 1; % save to graph
                    end
               
                else % horizontal
                    if d > gamma
                        graph(i,j) = 1; % save to graph
                    end
         
                end
            
            % start is not A, end is not B----------------------------------
            else
                if  data(j,5) == 1 % vertical
                    gamma = 1.5e4;
                    if d > gamma
                        graph(i,j) = 1; % save to graph
                    end
               
                else % horizontal             
                    gamma = 2e4;
                    if d > gamma
                        graph(i,j) = 1; % save to graph
                    end
         
                end
                
                               
            end % if end B
        end % if start A 
        
        end % end ix < jx
               
    end % end for
end % end for 

numnode = sum(graph,2)

%********************prepare data*****************************
function [data, datac] = data_prep(filename,flag)
    
data = csvread(filename);   
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
text(A(2),A(3),A(4),'  A');
text(B(2),B(3),B(4),'  B');
xlabel('x');
ylabel('y');
zlabel('z');
%-------------????---------------
%A???
xt = X(1);
yt = Y(1);
zt = Z(1);

x = X - xt;
y = Y - yt;
z = Z - zt;

position = [x y z];

% -----------????-------------
%???B???
xb = x(end);
yb = y(end);
zb = z(end);

cosPHI = xb/sqrt(xb^2+zb^2);
sinPHI = zb/sqrt(xb^2+zb^2);
cosTHETA = xb/sqrt(xb^2+yb^2);
sinTHETA = yb/sqrt(xb^2+yb^2);

%?z???THETA
Rz = [cosTHETA sinTHETA 0; -sinTHETA cosTHETA 0; 0 0 1]; %????Rz
%position1 = [cosTHETA sinTHETA 0; -sinTHETA cosTHETA 0; 0 0 1]* position';

%?y???PHI
Ry = [cosPHI 0 sinPHI; 0 1 0; -sinPHI 0 cosPHI]; %????Ry
%position2 = [cosPHI 0 sinPHI; 0 1 0; -sinPHI 0 cosPHI] * position1;

pos = Ry*Rz*position'; %??????????
pos(:,end); %?????B???

datac = [N pos' T L];
Ac = datac(1,:);
Bc = datac(end,:);
datacP = sortrows(datac(2:end-1,:),5);
figure;
scatter3(Ac(2),Ac(3),Ac(4),'r','o','filled'); %A
hold on;
scatter3(Bc(2),Bc(3),Bc(4),'r','o','filled'); %B
hold on;
scatter3(datacP(1:flag,2),datacP(1:flag,3),datacP(1:flag,4),'.','m');
hold on;
scatter3(datacP(flag+1:end,2),datacP(flag+1:end,3),datacP(flag+1:end,4),'.','b');
text(Ac(2),Ac(3),Ac(4),'  A');
text(Bc(2),Bc(3),Bc(4),'  B');
xlabel('x');
ylabel('y');
zlabel('z');

end