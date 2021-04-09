%% a
N=30;
teta_1=pi/10.25;
teta_2=2*pi/5;

n=0:1:N-1;
k_30=0:1:N-1;

s=2*cos(teta_1*n);
v=3*sin(teta_2*n);
x=s+v;

figure
plot(n,s)
hold on
plot(n,v)
hold on
plot(n,x)
title('Question 1 a time domain'); xlabel('n ','FontSize',12); legend([{'s[n]'};{'v[n]'};{'x[n]'}]);

X=fft(x);
V=fft(v);
S=fft(s);

figure
%stem(k,abs(X))
hold on
stem(k_30,abs(S))
hold on
stem(k_30,abs(V))
title('Question 1 a freq domain'); xlabel('k ','FontSize',12); legend([{'V[k]'};{'S[k]'}]);
%% b
k_45=0:1:44;
k_30=0:1:N-1;

x_z=[x zeros(1,45-length(x))];
X_z=fft(x_z);

figure
stem((2*pi/45).*k_45,abs(X_z),'LineWidth',1)
hold on
stem((2*pi/30).*k_30,abs(X),'LineStyle','--','LineWidth',1.5)
    
title('Question 1 b freq domain'); xlabel('normalized freq ','FontSize',12); legend([{'X_z[k]'};{'X[k]'}]);
%% c
n_2=0:1:44;

x_2=2*cos(teta_1*n_2)+3*sin(teta_2*n_2);
X_2=fft(x_2);

figure
stem(k_45,abs(X_2),'LineWidth',1);
hold on
stem(k_45,abs(X_z),'LineStyle','--','LineWidth',1);
    
title('Question 1 c freq domain'); xlabel('k ','FontSize',12); legend([{'X_2[k]'};{'X_z[k]'}]);

%% d
norm_X=norm(X)^2;
norm_x=norm(x)^2;
disp(norm_X/norm_x)

norm_X_z=norm(X_z)^2;
norm_x_z=norm(x_z)^2;
disp(norm_X_z/norm_x_z)
%% e

teta=0:0.01:2*pi;
n=[0 1 2 ];
h=[1/3 1/3 1/3];

H_f = exp(-1i*teta'*n)*h.' ; % dtft

figure
stem((2*pi/30).*k_30,abs(X)/30,'LineStyle','--','LineWidth',1)
hold on
plot(teta,abs(H_f))
title('Question 1 e freq domain'); xlabel('teta ','FontSize',12); legend([{'X[k]'};{'H_1(teta)'}]);


h_p=[h zeros(1,32-length(h))];
H_p=fft(h_p);
x_p=[x zeros(1,32-length(x))];
X_p=fft(x_p);
Y_1=X_p.*H_p;

figure
k_32=0:1:31;
stem((2*pi/32).*k_32,abs(Y_1),'LineStyle','--','LineWidth',1)
hold on
stem((2*pi/30).*k_30,abs(X),'LineWidth',1)
title('Question 1 e freq domain'); xlabel('normalized k ','FontSize',12); legend([{'Y_1[k]'};{'X[k]'}]);

y_1=ifft(Y_1);
y_1=y_1(1:30);

figure
plot(k_30,x,'LineWidth',1)
hold on
plot(k_30,s,'LineStyle','--')
hold on
plot(k_30,v,'LineStyle','--')
hold on
plot(k_30,y_1,'LineWidth',1)
title('Question 1 e time domain'); xlabel('n','FontSize',12); legend([{'x[t]'};{'s[t]'};{'v[k]'};{'y_1[k]'}]);

%% f
n=[0 1];
h_2=[1 1];

H_2_f = exp(-1i*teta'*n)*h_2.' ; % dtft

figure
stem((2*pi/30).*k_30,abs(X)/30,'LineStyle','--','LineWidth',1)
hold on
plot(teta,abs(H_2_f))
title('Question 1 f freq domain'); xlabel('teta ','FontSize',12); legend([{'X[k]'};{'H_2(teta)'}]);


h_2_p=[h_2 zeros(1,31-length(h_2))];
H_2_p=fft(h_2_p);
x_p_2=[x zeros(1,31-length(x))];
X_p_2=fft(x_p_2);
Y_2=X_p_2.*H_2_p;

figure
k_31=0:1:30;
stem((2*pi/31).*k_31,abs(Y_2),'LineStyle','--','LineWidth',1)
hold on
stem((2*pi/30).*k_30,abs(X),'LineWidth',1)
title('Question 1 f freq domain'); xlabel('normalized k ','FontSize',12); legend([{'Y_2[k]'};{'X[k]'}]);

y_2=ifft(Y_2);
y_2=y_2(1:30);

figure
plot(k_30,x,'LineWidth',1)
hold on
plot(k_30,s,'LineStyle','--')
hold on
plot(k_30,v,'LineStyle','--')
hold on
plot(k_30,y_2,'LineWidth',1)
title('Question 1 f time domain'); xlabel('n','FontSize',12); legend([{'x[t]'};{'s[t]'};{'v[k]'};{'y_2[k]'}]);

