
function eig_a1 = fun_prony(ya, dT, m)
%clear; clc; 
%load ymat

%dT = 0.001;  m=10;
% multiple-channel prony analysis
% y: each row a channel
% dT: sampling rate
% m:  system order
[N, n_ch] = size(ya);


% first, form Hankel matrix.
np = floor(N/3);
%np = 20; 
%m = 3; 
if(m > np)
   disp('order should be less than N/3');  
end

D=[]; Y=[];
for i=1:n_ch
[D1, Y1]=fun_prony_DY(ya(:,i),np); 
D =[D; D1];
Y =[Y; Y1];
end

[U,S,V]=svd(D);
D_prime = U(:,1:m)*S(1:m,1:m)*V(:,1:m)';
%D_prime = D;
a1= pinv(D_prime)*Y;

% Polynomial Roots Calculation
z = fun_a2R(a1);
eig_a1 = log(z)/dT; 

%show match by reconstructing signals. 
if(n_ch<= 5)
    row_plot = n_ch;
    col_plot = 1; 
else
    if(mod(sqrt(n_ch),1)>0)
        row_plot= floor(sqrt(n_ch))+1; 
    else
        row_plot=sqrt(n_ch);
    end
    if (mod(n_ch/row_plot,1)>0)
        col_plot = floor(n_ch/row_plot) + 1; 
    else 
        col_plot = n_ch/row_plot; 
    end
    % e.g., 6 signals: 3*2
    % e.g.; 7 signals: 3*3
end

%% signal reconstruction 
for i1=1:N;
    for j1=1:np;
        Z(i1,j1)=z(j1)^(i1-1);
    end 
end
for i=1:n_ch 
residue1(:,i) = pinv(Z)*ya(:,i);
y_hat(:,i)=Z*residue1(:,i);
end

figure('name', 'Prony Validation'); 
for i = 1: n_ch
subplot(row_plot, col_plot, i);  
%plot(dT*(1:N), y(:,i),'k', 'LineWidth',1);hold on;
plot(dT*(1:N), ya(:,i),'LineWidth',2);
hold on;
plot(dT*(1:N), real(y_hat(:,i)),'r' ); 
legend('Original', 'Prony');  grid on; 
end
 
figure('name', 'eigen plot');
scatter(real(eig_a1), imag(eig_a1)/2/pi,'Linewidth',2); 
ylabel('Hz'); xlabel('Real');
grid  on; 