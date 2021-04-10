function R = fun_a2R(a)
% given a vector a, find roots
[n,m] = size(a);
for k=1:n
   %s1(1,k)=a(n-k+1);
    s1(1,k) = a(k); 
end
p1 = [1 -s1];


% calculate the roots of plynomial 
R = roots(p1);
return
