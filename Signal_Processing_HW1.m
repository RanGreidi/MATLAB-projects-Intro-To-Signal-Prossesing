%% Question 1
%a
wm=3*pi;
t=0.2:0.01:3;

x = (4./(wm*pi*(t.^2))).*(sin(wm*t)).^2.*(cos(wm*t)).*(sin(2*wm*t));
figure;
plot(t,abs(x))
title(' Question 1a -Signal x(t) -time domain'); xlabel('t [sec]','FontSize',12); ylabel('x(t) [V]','FontSize',12);
hold off
%b


w= -17*pi:0.01:17*pi;
X = (2/j)*(triangularPulse((w-wm)./(4*wm))-triangularPulse((w+wm)./(4*wm)));
%X = (1/j)*(triangularPulse((w+wm)/(2*wm))-triangularPulse((w-wm)/(2*wm))+triangularPulse((w+3*wm)/(2*wm))-triangularPulse((w-3*wm)/(2*wm)));
figure;
plot(w,abs(X))
title('Question 1b -Signal X(w) -freq domain'); xlabel('w [rad/sec]','FontSize',12); ylabel('X(w) ','FontSize',12);
hold off
%c

Ts=1/20;
t_zoh=0.2:Ts:3;
ws=2*pi/Ts;

x_zoh=(4./(wm*pi*t_zoh.^2)).*(sin(wm*t_zoh)).^2.*(cos(wm*t_zoh)).*(sin(2*wm*t_zoh));

figure;
plot(t,x)
hold on ;    
stairs(t_zoh,x_zoh)
title(' Question 1c -Signal x(t) and x_Z_O_H(t) -time domain'); xlabel('t [sec]','FontSize',12); ylabel('x(t) [V]','FontSize',12);
legend([{'x(t)'};{'x_Z_O_H(t)'}]);
hold off

%d

X_p= zeros(1,length(w));
for k=-2:1:2
    X_p_temp=(1/Ts)*( (2/j)*(triangularPulse(((w-k*ws)-wm)./(4*wm))-triangularPulse(((w-k*ws)+wm)./(4*wm))) );
    X_p=X_p_temp +X_p;
end

X_zoh= X_p.*exp(-(j*(Ts/2)).*w).*sinc(w/(2*pi/Ts));

figure
%plot(w,abs(X_p))
plot(w,abs(X_zoh));


title(' Question 1d: Signal X_Z_O_H(w) -freq domain'); xlabel('w [rad/sec]','FontSize',12); ylabel('X_Z_O_H(w) ','FontSize',12);

%e
%X_rec=X_zoh.*(exp((j*(Ts/2)).*w))./(sinc(w/(2*pi/Ts))) ;

X_rec=X_zoh.*(exp((j*(Ts/2)).*w))./(sinc(w/(2*pi/Ts))).*rectangularPulse(w./ws) ;
x = (4./(wm*pi*(t.^2))).*(sin(wm*t)).^2.*(cos(wm*t)).*(sin(2*wm*t));
%X_rec= X_zoh ;

x_rect=[];
for i=0.2:0.01:3
    Y=(1/2*pi).*X_rec.*exp(j*w*i);
    ifur_for_para_t=trapz(w,Y);
    x_rect = [x_rect ifur_for_para_t];
end

x_rect=x_rect*(34*pi/(2*10682));  %normelize the trapz

%figure
%plot(w,abs(X_rec))


figure
plot(t,x_rect)
hold on
plot(t,x)
title(' Question 1e: Signal x_r_e_c(t)and X -time domain'); xlabel('t [sec]','FontSize',12);
legend([{'x_r_e_c(t)'};{'x(t)'}]);
hold off
%f
Ts=1/9;
ws=2*pi/Ts;

X_p= zeros(1,length(w));
for k=-2:1:2
    X_p_temp=(1/Ts)*( (2/j)*(triangularPulse(((w-k*ws)-wm)./(4*wm))-triangularPulse(((w-k*ws)+wm)./(4*wm))) );
    X_p=X_p_temp +X_p;
end
%figure
%plot(w,abs(X_p))

X_zoh= X_p.*exp(-(j*(Ts/2)).*w).*sinc(w/(2*pi/Ts));

X_rec=X_zoh.*(exp((j*(Ts/2)).*w))./(sinc(w/(2*pi/Ts))).*rectangularPulse(w./ws) ;

x_rect=[];
for i=0.2:0.01:3
    Y=(1/2*pi).*X_rec.*exp(j*w*i);
    ifur_for_para_t=trapz(w,Y);
    x_rect = [x_rect ifur_for_para_t];
end

x_rect=x_rect*(34*pi/(2*10682));  %normelize the trapz

figure
plot(w,abs(X_rec))
title(' Question 1f: Signal X_r_e_c(w)and x for Ts = 1/9 -freq domain'); xlabel('w [rad/sec]','FontSize',12);


figure
plot(t,x_rect)
hold on
plot(t,x)
title(' Question 1f: Signal x_r_e_c(t)and x for Ts = 1/9 -time domain'); xlabel('t [sec]','FontSize',12);
legend([{'x_r_e_c(t) Ts = 1/9'};{'x(t)'}]);
hold off


%% ------------------------------------------------------Question 2---------------------------------------------------------------------
%% Question 2
%a
t=0:0.01:2 ;
wa=7*pi ;
wb=5*pi ;
x=5*cos(wa*t)-3*sin(wb*t) ;
ts=0:2/15:2-2/15 ;  %fix it, it has only 14 samples the las one is the same as the first
w0=pi;
M=7;

x_semp=5*cos(wa*ts)-3*sin(wb*ts);

figure
plot(t,x)
hold on
stem(ts,x_semp)
title('Signal x(t) and Signal x_s_e_p(t) '); xlabel('t [sec]','FontSize',12); ylabel('[V]','FontSize',12);
legend([{'x(t)'};{'x_s_e_p(t)'}]);

%b

F = zeros(length(ts),length(ts));
for n = 1:length(ts)
   for k = 1:length(ts)       
       F(n,k)=exp(j*w0*ts(n)*(k-M-1));
   end
end

a=inv(F)*transpose(x_semp);


%c

x_rec=zeros(size(t));
for k = -M:1:M                                   % crrate the exponents for furrier series
   x_rec = x_rec + a(M+k+1).*exp(j*w0*k*t);
end

figure
plot(t,x_rec,'red','LineWidth',1.5)
hold on
plot(t,x,'--','LineWidth',1.5,'MarkerSize',10)
title('Signal x(t) and Signal x_r_e_c(t) '); xlabel('t [sec]','FontSize',12); ylabel('[V]','FontSize',12);
legend([{'x(t)'};{'x_r_e_c(t)'}]);

%d

ts_rand=2*rand(1,15);

x_semp_rand=5*cos(wa*ts_rand)-3*sin(wb*ts_rand);
%    repeat a
figure
plot(t,x)
hold on
stem(ts_rand,x_semp_rand)
title('Signal x(t) and Signal x_R_a_n_d_o_m_S_a_m_p(t)'); xlabel('t [sec]','FontSize',12); ylabel('[V]','FontSize',12);
legend([{'x(t)'};{'x_R_a_n_d_o_m_S_a_m_p(t)'}]);

%    repeat b
F = zeros(length(ts_rand),length(ts_rand));
for n = 1:length(ts_rand)
   for k = 1:length(ts_rand)       
       F(n,k)=exp(j*w0*ts_rand(n)*(k-M-1));
   end
end

a=inv(F)*transpose(x_semp_rand);

%    repeat c
x_rec_samples=zeros(size(t));
for k = -M:1:M                                                        % create the exponents for furrier series
   x_rec_samples = x_rec_samples + a(M+k+1).*exp(j*w0*k*t);
end

figure
plot(t,x_rec_samples,'red','LineWidth',1.5)
hold on
plot(t,x,'--','LineWidth',1.5,'MarkerSize',10)
title('Signal x(t) and Signal x_R_e_c_S_a_m_p_l_e_d(t) '); xlabel('t [sec]','FontSize',12); ylabel('[V]','FontSize',12);
legend([{'x(t)'};{'x_R_e_c_S_a_m_p_l_e_d(t)'}]);


%e

%-----repeat a, b, c

F = zeros(length(ts),length(ts));
for n = 1:length(ts)
   for k = 1:length(ts)       
       F(n,k)=exp(j*w0*(ts(n)+0.01*rand(1))*(k-M-1));
   end
end

a=inv(F)*transpose(x_semp);
%disp(a)

x_rec=zeros(size(t));
for k = -M:1:M                                   % creates the exponents for furrier series
   x_rec = x_rec + a(M+k+1).*exp(j*w0*k*t);
end

figure
plot(t,x_rec,'red','LineWidth',1.5)
hold on
plot(t,x,'--','LineWidth',1.5,'MarkerSize',10)
title('Signal x(t) and Signal x_r_e_c(t) but F with randon noise '); xlabel('t [sec]','FontSize',12); ylabel('[V]','FontSize',12);
legend([{'x(t)'};{'x_r_e_c(t)'}]);

cond_equal_samples = cond(F);

%-----repeat d

F = zeros(length(ts_rand),length(ts_rand));
for n = 1:length(ts_rand)
   for k = 1:length(ts_rand)       
       F(n,k)=exp(j*w0*(ts_rand(n)+0.01*rand(1))*(k-M-1));
   end
end

a=inv(F)*transpose(x_semp_rand);
disp(a)

x_rec_samples_RANDOM_F=zeros(size(t));
for k = -M:1:M                                                        % create the exponents for furrier series
   x_rec_samples_RANDOM_F = x_rec_samples_RANDOM_F + a(M+k+1).*exp(j*w0*k*t);
end

cond_random_samples = cond(F);

figure
plot(t,x,'red','LineWidth',1.5)
hold on
plot(t,x_rec_samples_RANDOM_F,'--','LineWidth',1.5,'MarkerSize',10)
title('Signal x(t) and Signal x_R_e_c_S_a_m_p_l_e_d but F with randon noise'); xlabel('t [sec]','FontSize',12); ylabel('[V]','FontSize',12);
legend([{'x(t)'};{'x_R_e_c_S_a_m_p_l_e_d(t)'}]);



% f

ts_rand_f=2*rand(1,40);
x_semp_rand_f=5*cos(wa*ts_rand_f)-3*sin(wb*ts_rand_f);

F = zeros(length(ts_rand_f),2*M+1);
for n = 1:length(ts_rand_f)
   for k = 1:2*M+1       
       F(n,k)=exp(j*w0*(ts_rand_f(n)+0.01*rand(1))*(k-M-1));
   end
end

F_f= pinv(F);

a_new=F_f*transpose(x_semp_rand_f);
disp(a_new)

x_rec_samples_RANDOM_F=zeros(size(t));
for k = -M:1:M                                                        % create the exponents for furrier series
   x_rec_samples_RANDOM_F = x_rec_samples_RANDOM_F + a_new(M+k+1).*exp(j*w0*k*t);
end

cond_random_samples_f = cond(F);

figure
plot(t,x,'red','LineWidth',1.5)
hold on
plot(t,x_rec_samples_RANDOM_F,'--','LineWidth',1.5,'MarkerSize',10)

title('Signal x(t) and Signal x_R_e_c_S_a_m_p_l_e_d but F with randon noise 40 samples'); xlabel('t [sec]','FontSize',12); ylabel('[V]','FontSize',12);
legend([{'x(t)'};{'x_R_e_c_S_a_m_p_l_e_d(t)'}]);


%% ----------------------------------------------Question 3-----------------------------------------------
%%
t=0:0.01:10;
T=10;

f=transpose(4*cos((4*pi/T)*t)+sin((10*pi/T)*t));
g=transpose(2*sign(sin((6*pi/T)*t))-4*sign(sin((4*pi/T)*t)));

n=-20:1:20;
phi_n=exp(j*(2*pi/T)*transpose(t).*n);   %makes matrix

ksi_k=[];                                %makes matrix
for k=0:1:19
    ksi_k = [ksi_k transpose(rectangularPulse((t/(T/20))-(k+0.5)))];
end

g_inr_pru_ksi=co(g,ksi_k,T);
s='coefficent for g with ksi set';
disp(s)
disp(g_inr_pru_ksi)

f_inr_pru_ksi=co(f,ksi_k,T);
s='coefficent for f with ksi set';
disp(s)
disp(f_inr_pru_ksi)


g_inr_pru_phi=co(g,phi_n,T);
s='coefficent for g with phi set';
disp(s)
disp(g_inr_pru_phi)

f_inr_pru_phi=co(f,phi_n,T);
s='coefficent for f with phi set';
disp(s)
disp(f_inr_pru_phi)


f_rec_ksi=zeros(size(t));  %reconsruct f by ksi set
for n=1:1:20
    f_rec_ksi=f_rec_ksi+f_inr_pru_ksi(n)*transpose(ksi_k(:,n));
end

g_rec_ksi=zeros(size(t));  %reconsruct g by ksi set
for n=1:1:20
    g_rec_ksi=g_rec_ksi+g_inr_pru_ksi(n)*transpose(ksi_k(:,n));
end

f_rec=zeros(size(t));  %reconsruct f by phi set
for n=1:1:41
    f_rec=f_rec+f_inr_pru_phi(n)*transpose(phi_n(:,n));
end

g_rec=zeros(size(t));  %reconsruct g by phi set
for n=1:1:41
    g_rec=g_rec+g_inr_pru_phi(n)*transpose(phi_n(:,n));
end



figure
plot(t,f,'red','LineWidth',1.5)
hold on
plot(t,f_rec,'--','LineWidth',1.5,'MarkerSize',10)
title('Signal f(t) and the signal f_r_e_c reconsructed by {phi} orthonormal set '); xlabel('t [sec]','FontSize',12); ylabel('[V]','FontSize',12);
legend([{'f(t)'};{'f_R_e_c(t)'}]);

figure
plot(t,g,'red','LineWidth',1.5)
hold on
plot(t,g_rec,'--','LineWidth',1.5,'MarkerSize',10)
title('Signal g(t) and the signal g_r_e_c reconsructed by {phi} orthonormal set '); xlabel('t [sec]','FontSize',12); ylabel('[V]','FontSize',12);
legend([{'g(t)'};{'g_R_e_c(t)'}]);


figure
plot(t,f,'red','LineWidth',1.5)
hold on
plot(t,f_rec_ksi,'--','LineWidth',1.5,'MarkerSize',10)
title('Signal f(t) and the signal f_r_e_c reconsructed by {ksi} orthonormal set '); xlabel('t [sec]','FontSize',12); ylabel('[V]','FontSize',12);
legend([{'f(t)'};{'f_R_e_c(t)'}]);

figure
plot(t,g,'red','LineWidth',1.5)
hold on
plot(t,g_rec_ksi,'--','LineWidth',1.5,'MarkerSize',10)
title('Signal g(t) and the signal g_r_e_c reconsructed by {ksi} orthonormal set '); xlabel('t [sec]','FontSize',12); ylabel('[V]','FontSize',12);
legend([{'g(t)'};{'g_R_e_c(t)'}]);

function [c]=co(a,M,T)
    size_a=size(a);
    size_M=size(M);
    c=zeros(size_M(2),1);
    t_integration=0:0.01:T;
    for i= 1:size(c)
       temp_culom = transpose(M(:,i));       
       numi=trapz(t_integration,transpose(a).*conj(temp_culom));
       denumi=trapz(t_integration,temp_culom.*conj(temp_culom));
       numi=(T/(2*size_a(1)))*numi;
       denumi=(T/(2*size_a(1)))*denumi;
       c(i)=numi./denumi;
    end   
end  







