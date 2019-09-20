clear;close all;

[data1, datac1] = data_prep('data1.csv',306);
[data2, datac2] = data_prep('data2.csv',167);
 
% save data
writematrix(data1, 'data1.csv');
writematrix(datac1, 'datac1.csv');
writematrix(data2, 'data2.csv');
writematrix(datac2, 'datac2.csv');


[graph1,W1] = build_graph(datac1,25,15,20,25,30,0.001);
[graph2,W2] = build_graph(datac2,20,10,15,20,20,0.001);

% save graph
writematrix(graph1,'graph1.csv');
writematrix(graph2,'graph2.csv');

[path1, distance1] = cal_shortestpath(graph1, W1)
[path2, distance2] = cal_shortestpath(graph2, W2)
% D1 = sqrt(datac1(304,2)^2 + datac1(304,3)^2 + datac1(304,4)^2) ...
% + sqrt((datac1(344,2)-datac1(304,2))^2 + (datac1(344,3)-datac1(304,3))^2 + (datac1(344,4)-datac1(304,4))^2) ...
% + sqrt((datac1(172,2)-datac1(344,2))^2 + (datac1(172,3)-datac1(344,3))^2 + (datac1(172,4)-datac1(344,4))^2) ...
% + sqrt((datac1(327,2)-datac1(172,2))^2 + (datac1(327,3)-datac1(172,3))^2 + (datac1(327,4)-datac1(172,4))^2);


function [path, distance] = cal_shortestpath(affinity, W)
% *********************calculate shortest path**********************************    
N = length(affinity);
s = [];
t = [];
w = [];

for i=1:N
    for j=1:N
        if affinity(i,j)==1
            s = [s,i];
            t = [t,j];
            w = [w,W(i,j)];
        end
    end
end

G = digraph(s,t,w);
% figure;
% p = plot(G,'EdgeLabel',G.Edges.Weight);

%Dijkstra:"positive" Bellman-Ford:"mixed"
[path, distance] = shortestpath(G,1,N,'Method', "positive");
% highlight(p, path,'EdgeColor','red')


end


function [graph,W] = build_graph(data,alpha1,alpha2,beta1,beta2,theta,delta)
%*********************build graph*****************************************

n = length(data);
graph = zeros(n);
W = Inf(n);

alpha = min(alpha1,alpha2);
beta = min(beta1,beta2);
gammav = min(alpha1,alpha2)/delta;
gammah = min(beta1,beta2)/delta;
gammaB = theta/delta;

for i = 1:n
    for j = 1:n       
        % Vertical：1  or Horizontal：0 
        
        % prerequiste: if ix < jx then calculate(j in front of i) 
        if data(i,2) < data(j,2) && data(i,2) >= 0 && data(j,2) >= 0
                  
        if data(i,5) == 10 %start A===============================================
            d = sqrt((data(j,2))^2+(data(j,3))^2+(data(j,4))^2);
            
            % start is A, end is B-----------------------------------------
            if data(j,5) == 100
                if d < gammaB
                    graph(i,j) = 1; % save to graph
                    W(i,j) = d;
                end
                
            % start is A, end is not B-------------------------------------
            else
                switch data(j,5) % j v or h
                    case 1  % vertical
                        if d < gammav
                            graph(i,j) = 1; % save to graph
                            W(i,j) = d;
                        end
               
                    case 0  % horizontal
                        if d < gammah
                            graph(i,j) = 1; % save to graph
                            W(i,j) = d;
                        end
                    otherwise
                end
            end
            
        else %start not A=========================================================
            
            if data(i,5) == 1 % i vertical
                dv = max(sqrt((data(j,2)-data(i,2))^2+(data(j,3)-data(i,3)+alpha)^2+(data(j,4)-data(i,4))^2),sqrt((data(j,2)-data(i,2))^2+(data(j,3)-data(i,3)-alpha)^2+(data(j,4)-data(i,4))^2)); 
            else              % i horizontal
                dh = max(sqrt((data(j,2)-data(i,2))^2+(data(j,3)-data(i,3))^2+(data(j,4)-data(i,4)+beta)^2),sqrt((data(j,2)-data(i,2))^2+(data(j,3)-data(i,3))^2+(data(j,4)-data(i,4)-beta)^2));       
            end
                     
            % start is not A, end is B--------------------------------------
            if data(j,5) == 100
                %gamma =3e4;
                switch data(i,5) % i v or h
                    case 1  % vertical
                        if dv < gammaB
                            graph(i,j) = 1; % save to graph
                            W(i,j) = dv;
                        end
               
                    case 0  % horizontal
                        if dh < gammaB
                            graph(i,j) = 1; % save to graph
                            W(i,j) = dh;
                        end
                    otherwise
                end
            
            % start is not A, end is not B----------------------------------
            else

                if  data(i,5)==1 && data(j,5) == 1 % i, j vertical
                    %gamma = 1.5e4;
                    if dv < gammav  %%%%%%%%%%%%%%%%%%%%%%%%%%
                        graph(i,j) = 1; % save to graph
                        W(i,j) = dv;
                    end           
                elseif data(i,5)==1 && data(j,5) == 0 % i vertical j horizontal             
                    %gamma = 2e4;
                    if dv < gammah
                        graph(i,j) = 1; % save to graph
                        W(i,j) = dv;
                    end
                elseif data(i,5)==0 && data(j,5) == 1 % i horizontal j vertical
                    if dh < gammav
                        graph(i,j) = 1; % save to graph
                        W(i,j) = dh;
                    end                  
                else % i,j horizontal
                    if dh < gammah
                        graph(i,j) = 1; % save to graph
                        W(i,j) = dh;
                    end
         
                end 
                
            end % if end B
        end % if start A 
        
        end % end ix < jx
               
    end % end for
end % end for 

%sum(graph,2); % connection num

end


function [data, datac] = data_prep(filename,flag)
%********************prepare data*****************************  

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


cosTHETA = xb/sqrt(xb^2+yb^2);
sinTHETA = yb/sqrt(xb^2+yb^2);

%?z???THETA
Rz = [cosTHETA sinTHETA 0; -sinTHETA cosTHETA 0; 0 0 1] %????Rz
%position1 = [cosTHETA sinTHETA 0; -sinTHETA cosTHETA 0; 0 0 1]* position';


a = Rz*[xb;yb;zb]; % first trans

cosPHI = a(1)/sqrt(a(1)^2+a(3)^2);
sinPHI = a(3)/sqrt(a(1)^2+a(3)^2);

%?y???PHI
Ry = [cosPHI 0 sinPHI; 0 1 0; -sinPHI 0 cosPHI] %????Ry
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