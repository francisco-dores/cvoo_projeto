%tema 50 - Ximas 1

clear
close all
clc

g = 9.81;

[cond_ini, max_deflec, inert, wing, deriv] = def_model();
w0 = cond_ini.aa0*cond_ini.u0;

%com a aprox das folhas
% a=[
% deriv.xu deriv.xw -w0 -g*cos(cond_ini.tt0);
% deriv.zu deriv.zw cond_ini.u0 -g*sin(cond_ini.tt0);
% deriv.mu+deriv.mwp*deriv.zu deriv.mw+deriv.mwp*deriv.zw deriv.mq+deriv.mwp*cond_ini.u0 -deriv.mwp*g*sin(cond_ini.tt0);
% 0 0 1 0];
% b=[deriv.xde deriv.xdsp;deriv.zde 0;deriv.mde+deriv.mwp*deriv.zde deriv.mdsp;0 0];

%considerando todas os coeficientes
a=[
deriv.xu deriv.xw -w0 -g*cos(cond_ini.tt0) 0;
deriv.zu/(1-deriv.zwp) deriv.zw/(1-deriv.zwp) (cond_ini.u0+deriv.zq)/(1-deriv.zwp) -g*sin(cond_ini.tt0)/(1-deriv.zwp) 0;
(deriv.mu+deriv.mwp*deriv.zu/(1-deriv.zwp)) (deriv.mw+deriv.mwp*deriv.zw/(1-deriv.zwp)) (deriv.mq+deriv.mwp*(cond_ini.u0+deriv.zq)/(1-deriv.zwp)) -deriv.mwp*g*sin(cond_ini.tt0)/(1-deriv.zwp) 0;
0 0 1 0 0;
0 -1 0 cond_ini.u0 0];

%u=[de;df;dsp]
b=[deriv.xde deriv.xdf deriv.xdsp;deriv.zde/(1-deriv.zwp) deriv.zdf/(1-deriv.zwp) deriv.zdsp/(1-deriv.zwp);deriv.mde+deriv.mwp*deriv.zde/(1-deriv.zwp) deriv.mdf+deriv.mwp*deriv.zdf/(1-deriv.zwp) deriv.mdsp+deriv.mwp*deriv.zdsp/(1-deriv.zwp);0 0 0;0 0 0];

c = eye(size(a));

d = zeros(size(b));

sys = ss(a,b,c,d);

[num_de,den_de]=ss2tf(a,b,c,d,1);

[num_df,den_df]=ss2tf(a,b,c,d,2);

[num_dsp,den_dsp]=ss2tf(a,b,c,d,3);

%sys = tf(num(1,:),den);
%rlocus(sys)

pzmap(sys)

damp(a)
[wn,zeta,p]=damp(a);

p_re = nonzeros(real(p));

%fugoide completo:
[M_fug,I_fug]=min(abs(p_re));
t_fug=log(2)/M_fug;

%período curto
[M_pc,I_pc]=max(abs(p_re));
t_pc=log(2)/M_pc;

disp('fugóide:')
if(p_re(I_fug)>0)
    disp(strcat('T_2 = ',num2str(t_fug)))
elseif(p_re(I_fug)<0)
    disp(strcat('T_{1/2} = ',num2str(t_fug)))
end

disp('período curto:')
if(p_re(I_pc)>0)
    disp(strcat('T_2 = ',num2str(t_pc)))
elseif(p_re(I_pc)<0)
    disp(strcat('T_{1/2} = ',num2str(t_pc)))
end

%% fugoide aproximado 
%(com diferentes variaveis pq o polinomio caracteristico
%de um sistema nao depende das variaveis de estado que o definem)

% a_fug=[
% deriv.xu -g*cond_ini.tt0;
% -deriv.zu/cond_ini.u0 0];

%-w0*g*sin(cond_ini.tt0)/(deriv.zq*cond_ini.u0)-g*cos(cond_ini.tt0)

a_fug=[
deriv.xu+deriv.zu*w0/(deriv.zq*cond_ini.u0) -w0*g*sin(cond_ini.tt0)/(deriv.zq*cond_ini.u0)-g*cos(cond_ini.tt0);
-deriv.zu/(deriv.zq*cond_ini.u0) g*sin(cond_ini.tt0)/(deriv.zq*cond_ini.u0)];

%b=[deriv.xde deriv.xdt;deriv.zde 0;deriv.mde+deriv.mwp*deriv.zde deriv.mdt;0 0];
%b_fug=[deriv.xde deriv.xdsp;-deriv.zde/cond_ini.u0 -deriv.zdsp/cond_ini.u0];

damp(a_fug) %fugoide estavel????


%% Período Curto Aproximado

a_pc=[
deriv.zw/(1-deriv.zwp) cond_ini.u0/(1-deriv.zwp);
deriv.mw+deriv.mwp*deriv.zw deriv.mq+deriv.mwp*cond_ini.u0];

damp(a_pc)


