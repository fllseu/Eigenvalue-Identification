clear;
clc;

load ISO-NE_case1

 figure (5) % Voltage magnitude
  plot(Time,Vm)
  xlabel('Time, s','FontSize',10)
  ylabel('KV','FontSize',10)
  title('Voltage magnitudes; phase-to-ground','FontSize',11)



y1=Vm(1202:1802,5);    %20 sec (from 40-60), 5 is number of the Vm signal
y2=Vm(1202:1802,13);

signal=[y1,y2]; 
Time1=Time(1202:1802); %20 sec



% figure (6)
% plot(Time1,[y1,y2])


dT=0.0344;  % dT should be 0.033-0.034


%%%%%%%%
%Set up the three methods
%%%%%%%%
%choose the mutiple signals and the system order
signal1=y1; 
signal2=y2;
order=3; %For all three methods
choose_L=200;  %For MP and ERA, choose L to be 200
            

               
choose_np_Prony=280; %For Prony, choose np to be 280

                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Prony
np=choose_np_Prony;
m=order;
ya1=[signal1,signal2];


%MP
M=order;
ya=[signal1,signal2];



%ERA
n=order;
h=[signal1,signal2].';
def=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rank Reduced Prony
%%%%%%%%%%%%%%%%%%%%%%%%%%
[N, n_ch] = size(ya1);
t = 0:dT:(N-1)*dT;


% first, form Hankel matrix.
%np = floor(N/3);
%np = 20; 
%m = 3; 
if(m > np)
   disp('order should be less than N/3');  
end

D=[]; Y=[];
for i=1:n_ch
[D1, Y1]=fun_prony_DY(ya1(:,i),np); 
D =[D; D1];
Y =[Y; Y1];
end

[U,S,V]=svd(D);
D_prime = U(:,1:m)*S(1:m,1:m)*V(:,1:m)';
%D_prime = D;
a1= pinv(D_prime)*Y;

% Polynomial Roots Calculation
z_Prony = fun_a2R(a1);
eig_a1 = log(z_Prony)/dT; 

% reserve only m eigenvalues. 
[re, ind]=sort(real(eig_a1),'descend');

eig_s = eig_a1(ind(1:m));
z_Prony = exp(eig_a1(ind(1:m))*dT);


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
    for j1=1:m;
        Z_Prony(i1,j1)=z_Prony(j1)^(i1-1);
    end 
end
for i=1:n_ch 
residue1_Prony(:,i) = pinv(Z_Prony)*ya1(:,i);
y_hat_Prony(:,i)=Z_Prony*residue1_Prony(:,i);
end

   

% figure(101);
% label = {'Vm 2 (kV)','Vm 4 (kV)'}; %for multiple channel
% for i = 1: n_ch
% subplot(row_plot, col_plot, i);  
% %plot(dT*(1:N), y(:,i),'k', 'LineWidth',1);hold on;
% plot(t, ya1(:,i),'b','LineWidth',2);
% hold on;
% plot(t, real(y_hat_Prony(:,i)),'r--' ,'LineWidth',1); 
% grid on;
%     legend('Original', 'Prony');
%     ylabel(label(1,i)); %for multiple channel
% end
%  xlabel('Time (s)')




% figure(102);
% plot(real(eig_s)/2/pi, ((imag(eig_s)))/2/pi,'r+','Linewidth',2, 'Markersize',10);
% grid on;
% hold on;
% grid on;
% ax = gca;
% ax.GridLineStyle = ':';
% ax.GridAlpha = 0.5;
% title('Eigenvalues (Hz)')
% xlabel('Real')
% ylabel('Imaginary')
% % xlim([-37.77 -37.66]);
% legend('Original', 'Prony');


% % %%%%%%%%%%%%%%%%%%%%%%%%%
% % %MP
% % %%%%%%%%%%%%%%%%%%%%%%%%%
n_ch_MP = size(ya, 2); 
N_MP = size(ya,1)-1; t1 = 0:dT:N_MP*dT;

if(n_ch_MP<= 5)
    row_plot_MP = n_ch_MP;
    col_plot_MP = 1; 
else
    if(mod(sqrt(n_ch_MP),1)>0)
        row_plot_MP= floor(sqrt(n_ch_MP))+1; 
    else
        row_plot_MP=sqrt(n_ch_MP);
    end
    if (mod(n_ch_MP/row_plot_MP,1)>0)
        col_plot_MP = floor(n_ch_MP/row_plot_MP) + 1; 
    else 
        col_plot_MP = n_ch_MP/row_plot_MP; 
    end
    % e.g., 6 signals: 3*2
    % e.g.; 7 signals: 3*3
end

% figure(999);
% for i=1:n_ch_MP
%     subplot(row_plot_MP, col_plot_MP, i); plot(t1, ya(:,i),'b','linewidth',2); %legend('original signal');
%     hold on;
% end

D =[];
%N = size(ya,1)-1; t1 = 0:dT:N*dT;
% L = floor(1/3*N_MP);
L=choose_L;
%M=10; % order 

for k=1:size(ya,2) % visit each column
    % for each channel, build a Hankel matrix
    for i=1:L+2  
        H(i,:) = ya(i:i+N_MP-L-1); 
    end
    D =[D, H];
end

[U,S,V] = svd(D); 
%figure(999)
%semilogy(diag(S));

U_prime= U(:,1:M);
U1 = U_prime(1:L+1, :);
U2 = U_prime(2:end, :); 

Lambda_MP = inv(U2'*U1)*(U2'*U2);
z_MP=eig(Lambda_MP);
eig_s_MP = log(z_MP)/dT;

% figure(888);
% semilogy(real(eig_s_MP), (abs(imag(eig_s_MP))/2/pi),'+','Linewidth',2, 'Markersize',10);
% ylabel('Hz')
% title('eigenvalues'); grid on;
% figure(889);
% plot((real(eig_s_MP)/2/pi), ((imag(eig_s_MP))/2/pi),'+','Linewidth',2, 'Markersize',10);
% ylabel('Imaginary')
% title('eigenvalues (Hz)'); grid on;



%% signal reconstruction 
for i1=1:N_MP+1;
    for j1=1:M;
        Z_MP(i1,j1)=z_MP(j1)^(i1-1);
    end 
end
for i=1:size(ya,2) % for three signal, reconstruct
residue1_MP(:,i) = pinv(Z_MP)*ya(:,i);
y_hat_MP(:,i)=Z_MP*residue1_MP(:,i);
end

% figure(999); 
% label = {'Vm 2 (kV)','Vm 4 (kV)'}; %for multiple channel
% for i=1:n_ch_MP
%     subplot(row_plot_MP, col_plot_MP, i);      
% %     plot(t1, ya(:,i),'Linewidth',2); hold on;
%     plot(t1, real(y_hat_MP(:,i)),'r','Linewidth',1);  
%     grid on
%     legend('Original', 'MP');
%         ylabel(label(1,i)); %for multiple channel
% end
% xlabel('Time (s)')
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %ERA
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
[n_ch_ERA, N1_ERA] = size(h);
N_ERA = N1_ERA-2; 
% L_ERA = floor(N_ERA/3); 
L_ERA = choose_L; 
H0=[]; H1=[];
for k = 1: N_ERA+1-L_ERA
    H0 = [H0; h(:,k+1:k+L_ERA)];
    H1 = [H1; h(:,k+2: k+1+L_ERA)]; 
end


% Factorization of the Hankel matrix by use of SVD
[R_ERA,Sigma_ERA,S_ERA] = svd(H0);   
[U1_ERA,S1_ERA,V1_ERA] = svd(H1);   
H1_v1 = U1_ERA(:,1:n)*S1_ERA(1:n, 1:n)*V1_ERA(:,1:n)'; 
% R and S are orthonormal and Sigma is a rectangular matrix

Sigman = Sigma_ERA(1:n,1:n);            

Wo = R_ERA(:,1:n)*Sigman^0.5;           % observability matrix
Co = Sigman^.5*S_ERA(:,1:n)';           % controllability matrix

% The identified system matrix are:
A = Sigman^-.5*R_ERA(:,1:n)'*H1_v1*S_ERA(:,1:n)*Sigman^-.5;            % dynamic matrix
B = Co(:,1);                    % input matrix
C = Wo(1,:);                    % output matrix
D = h(1);                       % direct-transmission matrix

sysdisc = ss(A,B,C,D,dT);       % discrete-time system

if def == 2                            
    syscont = d2c(sysdisc,'zoh');       % Conversion of discrete LTI models to continuous time
    [A,B,C,D]=ssdata(syscont);          % continuous system
end

%--------------------------------------------------------------------------
% show signal match
% show eigenvalue
z_ERA = eig(A); 
eig_s_ERA = log(eig(A))/dT; 
% figure(888); 
% scatter((real(eig_s_ERA)/2/pi),(imag(eig_s_ERA)/2/pi),'LineWidth',2); 
% title('eigenvalues (Hz)');
% grid on; 

%show match by reconstructing signals. 
if(n_ch_ERA<= 5)
    row_plot_ERA = n_ch_ERA;
    col_plot_ERA = 1; 
else
    if(mod(sqrt(n_ch_ERA),1)>0)
        row_plot_ERA= floor(sqrt(n_ch_ERA))+1; 
    else
        row_plot_ERA=sqrt(n_ch_ERA);
    end
    if (mod(n_ch_ERA/row_plot_ERA,1)>0)
        col_plot_ERA = floor(n_ch_ERA/row_plot_ERA) + 1; 
    else 
        col_plot_ERA = n_ch_ERA/row_plot_ERA; 
    end
    % e.g., 6 signals: 3*2
    % e.g.; 7 signals: 3*3
end
% figure(999);

ya_ERA =h.'; 
%% signal reconstruction 
for i1=1:N_ERA+2;
    for j1=1:n;
        Z_ERA(i1,j1)=z_ERA(j1)^(i1-1);
    end 
end
for i=1:size(ya_ERA,2) % for three signal, reconstruct
residue1_ERA(:,i) = pinv(Z_ERA)*ya_ERA(:,i);
y_hat_ERA(:,i)=Z_ERA*residue1_ERA(:,i);
end

% figure(998); 
% label = {'Vm 2 (kV)','Vm 4 (kV)'}; %for multiple channel
% t1 = 0:dT:(N_ERA+1)*dT;
% for i=1:n_ch_ERA
%     subplot(row_plot_ERA, col_plot_ERA, i);      
%     plot(t1, ya_ERA(:,i),'b','Linewidth',2); hold on;
%     plot(t1, real(y_hat_ERA(:,i)),'r-.','Linewidth',1);  
%     grid on
%     legend('Original', 'ERA');
%     ylabel(label(1,i)); %for multiple channel
% end
% xlabel('Time (s)')
% 
% 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% 
figure(1);
subplot(2,1,1)
plot(t1, ya(:,1),'b','Linewidth',2);
hold on;
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.5;
plot(t1, real(y_hat_Prony(:,1)),'r--', 'linewidth', 1.4); 
hold on;
plot(t1, real(y_hat_MP(:,1)),'c','Linewidth',1.4); 
hold on;
plot(t1, real(y_hat_ERA(:,1)),'k-.','Linewidth',1.4); 
ylabel('Vm 5 (kV)')
title('Reconstructed signals against the original measurements')
subplot(2,1,2)
plot(t1, ya(:,2),'b','Linewidth',2);
hold on;
grid on;
plot(t1, real(y_hat_Prony(:,2)),'r--', 'linewidth', 1.4); 
hold on;
plot(t1, real(y_hat_MP(:,2)),'c','Linewidth',1.4); 
hold on;
plot(t1, real(y_hat_ERA(:,2)),'k-.','Linewidth',1.4);
ylabel('Vm 13 (kV)')
xlabel('Time (s)')
legend('Original', 'Prony', 'MP', 'ERA');


% % eig_orig=pole(sys2);
% % eigenvalues=[eig_orig, eig_a1, eig_s_MP, eig_s_ERA];
% % figure(2);
% % color_str ={'bd','r+','g*','ko'};
% % for i=1:size(eigenvalues,2);
% % plot(real(eigenvalues(:,i)), ((imag(eigenvalues(:,i)))),color_str{i},'Linewidth',2, 'Markersize',10);
% % grid on;
% % hold on;
% % grid on;
% % ax = gca;
% % ax.GridLineStyle = ':';
% % ax.GridAlpha = 0.5;
% % title('Eigenvalues (rad/s)')
% % xlabel('Real')
% % ylabel('Imaginary')
% % % xlim([-37.77 -37.66]);
% % legend('Original', 'Prony', 'MP', 'ERA');
% % end
% 
% 
% 
% 
% eig_orig=pole(sys2);
% figure(2);
% plot(real(eig_orig), imag(eig_orig),'bd','Linewidth',2, 'Markersize',10);
% hold on;
eigenvalues=[eig_s, eig_s_MP, eig_s_ERA];
figure(2);
% color_str ={'r+','g*','ko'};
color_str ={'r+','b*','ko'};
for i=1:size(eigenvalues,2);
plot(real(eigenvalues(:,i)), ((imag(eigenvalues(:,i)))),color_str{i},'Linewidth',2, 'Markersize',10);
grid on;
hold on;
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.5;
title('Eigenvalues (rad/s)')
xlabel('Real')
ylabel('Imaginary')
% xlim([-37.77 -37.66]);
legend('Prony', 'MP', 'ERA');
end


% figure(3);
% plot(real(eig_orig)/2/pi, imag(eig_orig)/2/pi,'bd','Linewidth',2, 'Markersize',10);
% hold on;
eigenvalues=[eig_s, eig_s_MP, eig_s_ERA];
figure(3);
% color_str ={'r*','g+','ko'};
color_str ={'r+','b*','ko'};

for i=1:size(eigenvalues,2);
plot(real(eigenvalues(:,i)/2/pi), ((imag(eigenvalues(:,i)))/2/pi),color_str{i},'Linewidth',2, 'Markersize',10);
grid on;
hold on;
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.5;
title('Eigenvalues (Hz)')
xlabel('Real')
ylabel('Imaginary')
% xlim([-37.77 -37.66]);
legend( 'Prony', 'MP', 'ERA');
end

% 
% figure(44);
% plot(real(eig_orig)/2/pi, imag(eig_orig)/2/pi,'bd','Linewidth',2, 'Markersize',10);
% hold on;
% eigenvalues=[eig_s_ERA];
% figure(44);
% color_str ={'ro'};
% for i=1:size(eigenvalues,2);
% plot(real(eigenvalues(:,i)/2/pi), ((imag(eigenvalues(:,i)))/2/pi),color_str{i},'Linewidth',2, 'Markersize',10);
% grid on;
% hold on;
% grid on;
% ax = gca;
% ax.GridLineStyle = ':';
% ax.GridAlpha = 0.5;
% title('Eigenvalues (Hz)')
% xlabel('Real')
% ylabel('Imaginary')
% % xlim([-37.77 -37.66]);
% legend('Original', 'ERA');
% end
% 
% 
% % figure(2);
% % grid on;
% % hold on;
% % plot(real(eig_s_MP), ((imag(eig_s_MP))),'r+','Linewidth',2, 'Markersize',10);
% % hold on;
% % plot(real(eig_s_MP), ((imag(eig_s_ERA))),'r+','Linewidth',2, 'Markersize',10);
% 
%  
% 


Prony_modes=imag(eig_s/2/pi);
ERA_modes=imag(eig_s_ERA/2/pi);
MP_modes=imag(eig_s_MP/2/pi);