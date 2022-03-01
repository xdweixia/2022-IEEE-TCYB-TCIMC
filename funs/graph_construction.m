function S=graph_construction(X,k)
v= size(X,2);
if nargin < 2
    for i=1:v
      num(i)=  size(X{i},2);
%         k(i) =num(i)-2;
        k(i) =9;
    end
end;
%initialize weighted_distX
for i = 1:v 
    distX_initial{i} =  L2_distance_1( X{i},X{i}) ;                  %initialize X
    [distXs, idx] = sort(distX_initial{i},2);
    S{i}=zeros(num(i),num(i));
    for j = 1:num(i)
        di = distXs(j,2:k(i)+2);
        id = idx(j,2:k(i)+2);
        S{i}(j,id) = (di(k(i)+1)-di)/(k(i)*di(k(i)+1)-sum(di(1:k(i)))+eps);               %initialize S
    end
end
for i=1:v
S{i}=(S{i}'+S{i})*0.5;
end