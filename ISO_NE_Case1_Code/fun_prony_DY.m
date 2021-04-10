function [D1,Y1]=fun_prony_DY(y,np)
% np is the order of the model
% y is the input

N1 = length(y); 
% n is the system order, N as data, Da = Y, Y is N, D dimension: N*n
% N1 = N+n. 
N = N1-np; 
n1 = np;
%form the first matrix [n*1]  'Y1'
Y1 = zeros(N,1);
for  k=1:N
     Y1(k,1)=y(n1+k);
end

%Form the D matrix
D1 = zeros(N,n1);

for i=1:N
    for j = 1:n1
       %D(i,j) = pron( (i-1)+(n-1)-(j-1));  % based on eq(2) in 2012 PWRS
       %paper Zhou, N, " A stepwise regression method for estimating dominant
       %Electromechanical modes"
       D1(i,j) = y((i-1)+(n1-1)-(j-1) +1); % start from y(1) instead of y(0)
    end
end
return
end

