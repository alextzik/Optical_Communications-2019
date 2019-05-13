%Constants
V=1.5*10^(-16);
N0=1.5*10^8;
gamma=0.4;
beta=0.01;
t_c=1*10^(-9);
t_p=3*10^(-12);
e_NL=3.33*10^(-7);
s0=2.5*10^(-20);
L=300*10^(-6);
c=3*10^8;
vg=c/3.527;
R1=0.3;
R2=0.3;
g0=vg*s0;
lambda0=1550*10^(-9);
e=1.60217662 * 10^(-19);
amir=log(1/(R1*R2))/(2*L);
h=6.62607004 * 10^(-34);
f=c/lambda0;
%%
%----------a. (DC characteristic equation of laser diode)-----------
Nth=V/(gamma*g0*t_p)+N0;
Ith=e*Nth/t_c;
% i)
I=0:0.001:5*Ith;
Ps_i=zeros(1, length(I));
i=1;
for I=0:0.001:5*Ith;
  if(I>Ith)
    Ps_i(i)=gamma*t_p/e*(I-Ith);
  end
  i=i+1;
end

Pe_i=0.5*vg*amir*h*f*Ps_i;
% ii)
syms N P

Ps_ii=zeros(1, length(0:0.001:5*Ith));
Ns_ii=zeros(1, length(0:0.001:5*Ith));
i=1;
for I=0:0.001:5*Ith
   [solN, solP]=solve(I/e-N/t_c-g0*(N-N0)/V*P==0, gamma*g0*(N-N0)/V*P-P/t_p==0);
   if (I<=Ith)
       Ps_ii(1,i)=subs(solP(find(subs(solP)==0)));
       Ns_ii(1,i)=subs(solN(find(subs(solP)==0)));
   else 
       Ps_ii(1,i)=subs(solP(find(subs(solP)>0)));
       Ns_ii(1,i)=subs(solN(find(subs(solP)>0)));
   end
   
   i=i+1;
end
Pe_ii=0.5*vg*amir*h*f*Ps_ii;

I=0:0.001:5*Ith;

%plots
figure(1);
plot(I, Pe_i);

figure(2);
plot(I, Pe_ii);

figure(3);
plot(I, Pe_i);
hold on;
plot(I, Pe_ii);

%%
%--------------b. (response to single pulse) -------------
e_NL=3.33*10^(-8);
t0=2*10^(-9);
p=@(t)(rectpuls(t-t0/2, t0));


%i)1)
I=@(t)(0+2*Ith*p(t));
[solN, solP]=solve(0/e-N/t_c-g0*(N-N0)/(V*(1+e_NL*P))*P==0, gamma*g0*(N-N0)/(V*(1+e_NL*P))*P-P/t_p+gamma*beta*N/t_c==0);
p_t=subs(solP(find(subs(solP)==0)));
n=subs(solN(find(subs(solP)==0)));
numPhot=double(p_t);
numCar=double(n);
%N=y(1), P=y(2)
dydt=@(t, y)[I(t)/e-y(1)/t_c-g0/V*(y(1)-N0)*y(2)/(1+e_NL*y(2)); gamma*g0*(y(1)-N0)/V*y(2)/(1+e_NL*y(2))-y(2)/t_p+gamma*beta*y(1)/t_c];
[t, y]=ode45(dydt, [-t0/2, 3*t0], [Ns_ii(1,1); Ps_ii(1,1)]);

figure(4);
Is=I(-t0/2:t0/100:3*t0);
hold on;
ax=plotyy(-t0/2:t0/100:3*t0, Is, t, 0.5*vg*amir*h*f*y(:,2));
xlabel('time (sec)');
ylabel(ax(1), 'I (current)');
ylabel(ax(2), 'Εξωτερικά Εκπεμπόμενη Ισχύς');
title('Εξωτερικά Εκπεμπόμενη Ισχύς');

%i)2)
%Find initial conditions from dc bias current
I=@(t)(Ith+Ith*p(t));
[solN, solP]=solve(Ith/e-N/t_c-g0*(N-N0)/(V*(1+e_NL*P))*P==0, gamma*g0*(N-N0)/(V*(1+e_NL*P))*P-P/t_p+gamma*beta*N/t_c==0);
p_t=subs(solP(find(subs(solP)>=0)));
n=subs(solN(find(subs(solP)>=0)));
numPhot=double(p_t);
numCar=double(n);

%N=y(1), P=y(2)
dydt=@(t, y)[I(t)/e-y(1)/t_c-g0/V*(y(1)-N0)*y(2)/(1+e_NL*y(2)); gamma*g0*(y(1)-N0)/V*y(2)/(1+e_NL*y(2))-y(2)/t_p+gamma*beta*y(1)/t_c];
[t, y]=ode45(dydt, [-t0/2, 3*t0], [numCar ; numPhot]);
output=0.5*vg*amir*h*f*y(:,2);

figure(5);
Is=I(-t0/2:t0/100:3*t0);
hold on;
ax=plotyy(-t0/2:t0/100:3*t0, Is, t, output);
xlabel('time (sec)');
ylabel(ax(1), 'I (current)');
ylabel(ax(2), 'Εξωτερικά Εκπεμπόμενη Ισχύς');
title('Εξωτερικά Εκπεμπόμενη Ισχύς');
%%
%ii)1)
beta=0;
%Find initial conditions from dc bias current
I=@(t)(1.1*Ith+1.4*Ith*p(t));
[solN, solP]=solve(1.1*Ith/e-N/t_c-g0*(N-N0)/(V*(1+e_NL*P))*P==0, gamma*g0*(N-N0)/(V*(1+e_NL*P))*P-P/t_p==0);
p_t=subs(solP(find(subs(solP)>0)));
n=subs(solN(find(subs(solP)>0)));
numPhot=double(p_t);
numCar=double(n);

%N=y(1), P=y(2)
dydt=@(t, y)[I(t)/e-y(1)/t_c-g0/V*(y(1)-N0)*y(2)/(1+e_NL*y(2)); gamma*g0*(y(1)-N0)/V*y(2)/(1+e_NL*y(2))-y(2)/t_p+gamma*beta*y(1)/t_c];
[t, y]=ode45(dydt, [-t0/2, 3*t0], [numCar ; numPhot]);

figure(6);
Is=I(-t0/2:t0/100:3*t0);
hold on;
ax=plotyy(-t0/2:t0/100:3*t0, Is, t, 0.5*vg*amir*h*f*y(:,2));
xlabel('time (sec)');
ylabel(ax(1), 'I (current)');
ylabel(ax(2), 'Εξωτερικά Εκπεμπόμενη Ισχύς');
title('Εξωτερικά Εκπεμπόμενη Ισχύς');

%ii)2)
beta=0;
e_NL=3.33*10^(-7);
%Find initial conditions from dc bias current
I=@(t)(1.1*Ith+1.4*Ith*p(t));
[solN, solP]=solve(1.1*Ith/e-N/t_c-g0*(N-N0)/(V*(1+e_NL*P))*P==0, gamma*g0*(N-N0)/(V*(1+e_NL*P))*P-P/t_p==0);
p_t=subs(solP(find(subs(solP)>0)));
n=subs(solN(find(subs(solP)>0)));
numPhot=double(p_t);
numCar=double(n);

%N=y(1), P=y(2)
dydt=@(t, y)[I(t)/e-y(1)/t_c-g0/V*(y(1)-N0)*y(2)/(1+e_NL*y(2)); gamma*g0*(y(1)-N0)/V*y(2)/(1+e_NL*y(2))-y(2)/t_p+gamma*beta*y(1)/t_c];
[t, y]=ode45(dydt, [-t0/2, 3*t0], [numCar ; numPhot]);

figure(7);
Is=I(-t0/2:t0/100:3*t0);
hold on;
ax=plotyy(-t0/2:t0/100:3*t0, Is, t, 0.5*vg*amir*h*f*y(:,2));
xlabel('time (sec)');
ylabel(ax(1), 'I (current)');
ylabel(ax(2), 'Εξωτερικά Εκπεμπόμενη Ισχύς');
title('Εξωτερικά Εκπεμπόμενη Ισχύς');

%%
%-----------------c. (Response to PRBS sequence) ----------------------------
beta=0;
e_NL=3.33*10^(-7);
prbs=[1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0];
nofPRBS=length(prbs);
%i)
Bt=2.5*10^9;
pulseDuration=1/Bt;
t = 0:pulseDuration/50:nofPRBS*pulseDuration;
d = [0:pulseDuration:(nofPRBS-1)*pulseDuration; prbs]';
x = @(var)rectpuls(var-pulseDuration/2, pulseDuration);
p = @(t) pulstran(t,d,x);
%Find initial conditions from dc bias current
I=@(t)(1.1*Ith+1.4*Ith*p(t));
[solN, solP]=solve(1.1*Ith/e-N/t_c-g0*(N-N0)/(V*(1+e_NL*P))*P==0, gamma*g0*(N-N0)/(V*(1+e_NL*P))*P-P/t_p==0);
p_t=subs(solP(find(subs(solP)>0)));
n=subs(solN(find(subs(solP)>0)));
numPhot=double(p_t);
numCar=double(n);

%N=y(1), P=y(2)
dydt=@(t, y)[I(t)/e-y(1)/t_c-g0/V*(y(1)-N0)*y(2)/(1+e_NL*y(2)); gamma*g0*(y(1)-N0)/V*y(2)/(1+e_NL*y(2))-y(2)/t_p+gamma*beta*y(1)/t_c];
[t, y]=ode45(dydt, [0, nofPRBS*pulseDuration], [numCar ; numPhot]);

%Pe
figure(8);
Is=I(t);
hold on;
ax=plotyy(t, Is, t, 0.5*vg*amir*h*f*y(:,2));
xlabel('time (sec)');
ylabel(ax(1), 'I (current)');
ylabel(ax(2), 'Εξωτερικά Εκπεμπόμενη Ισχύς');
title('Εξωτερικά Εκπεμπόμενη Ισχύς');

%N
figure(9);
ax=plotyy(t, Is, t, y(:,1)/Nth);
xlabel('time (sec)');
ylabel(ax(1), 'I (current)');
ylabel(ax(2), 'Κανονικοποιημένος Αριθμός Φορέων');
title('Κανονικοποιημένος Αριθμός Φορέων');
%Eye
figure(10);
j=1;
hold on;
for i=1:1:31
   Pe_eye=[];
   while(t(j)<pulseDuration*i)
         Pe_eye=[Pe_eye, 0.5*vg*amir*h*f*y(j,2)];
         j=j+1;
   end
   plot(0:pulseDuration/length(Pe_eye):pulseDuration-pulseDuration/length(Pe_eye), Pe_eye);
end
xlabel('time (sec)');
ylabel('Εξωτερικά Εκπεμπόμενη Ισχύς');
title('Διάγραμμα Οφθαλμού Εκπεμπόνενης Ισχύος');
%%
%ii)
Bt=5*10^9;
pulseDuration=1/Bt;
t = 0:pulseDuration/50:nofPRBS*pulseDuration;
d = [0:pulseDuration:(nofPRBS-1)*pulseDuration; prbs]';
x = @(var)rectpuls(var-pulseDuration/2, pulseDuration);
p = @(t) pulstran(t,d,x);
%Find initial conditions from dc bias current
I=@(t)(1.1*Ith+1.4*Ith*p(t));
[solN, solP]=solve(1.1*Ith/e-N/t_c-g0*(N-N0)/(V*(1+e_NL*P))*P==0, gamma*g0*(N-N0)/(V*(1+e_NL*P))*P-P/t_p==0);
p_t=subs(solP(find(subs(solP)>0)));
n=subs(solN(find(subs(solP)>0)));
numPhot=double(p_t);
numCar=double(n);

%N=y(1), P=y(2)
dydt=@(t, y)[I(t)/e-y(1)/t_c-g0/V*(y(1)-N0)*y(2)/(1+e_NL*y(2)); gamma*g0*(y(1)-N0)/V*y(2)/(1+e_NL*y(2))-y(2)/t_p+gamma*beta*y(1)/t_c];
[t, y]=ode45(dydt, [0, nofPRBS*pulseDuration], [numCar ; numPhot]);

%Pe
figure(11);
Is=I(t);
hold on;
ax=plotyy(t, Is, t, 0.5*vg*amir*h*f*y(:,2));
xlabel('time (sec)');
ylabel(ax(1), 'I (current)');
ylabel(ax(2), 'Εξωτερικά Εκπεμπόμενη Ισχύς');
title('Εξωτερικά Εκπεμπόμενη Ισχύς');
%N
figure(12);
ax=plotyy(t, Is, t, y(:,1)/Nth);
xlabel('time (sec)');
ylabel(ax(1), 'I (current)');
ylabel(ax(2), 'Κανονικοποιημένος Αριθμός Φορέων');
title('Κανονικοποιημένος Αριθμός Φορέων');
%Eye
figure(13);
j=1;
hold on;
for i=1:1:31
   Pe_eye=[];
   while(t(j)<pulseDuration*i)
         Pe_eye=[Pe_eye, 0.5*vg*amir*h*f*y(j,2)];
         j=j+1;
   end
   plot(0:pulseDuration/length(Pe_eye):pulseDuration-pulseDuration/length(Pe_eye), Pe_eye);
end
xlabel('time (sec)');
ylabel('Εξωτερικά Εκπεμπόμενη Ισχύς');
title('Διάγραμμα Οφθαλμού Εκπεμπόνενης Ισχύος');
%%

%iii)
Bt=10*10^9;
pulseDuration=1/Bt;
t = 0:pulseDuration/50:nofPRBS*pulseDuration;
d = [0:pulseDuration:(nofPRBS-1)*pulseDuration; prbs]';
x = @(var)rectpuls(var-pulseDuration/2, pulseDuration);
p = @(t) pulstran(t,d,x);
%Find initial conditions from dc bias current
I=@(t)(1.1*Ith+1.4*Ith*p(t));
[solN, solP]=solve(1.1*Ith/e-N/t_c-g0*(N-N0)/(V*(1+e_NL*P))*P==0, gamma*g0*(N-N0)/(V*(1+e_NL*P))*P-P/t_p==0);
p_t=subs(solP(find(subs(solP)>0)));
n=subs(solN(find(subs(solP)>0)));
numPhot=double(p_t);
numCar=double(n);

%N=y(1), P=y(2)
dydt=@(t, y)[I(t)/e-y(1)/t_c-g0/V*(y(1)-N0)*y(2)/(1+e_NL*y(2)); gamma*g0*(y(1)-N0)/V*y(2)/(1+e_NL*y(2))-y(2)/t_p+gamma*beta*y(1)/t_c];
[t, y]=ode45(dydt, [0, nofPRBS*pulseDuration], [numCar ; numPhot]);

%Pe
figure(14);
Is=I(t);
hold on;
ax=plotyy(t, Is, t, 0.5*vg*amir*h*f*y(:,2));
xlabel('time (sec)');
ylabel(ax(1), 'I (current)');
ylabel(ax(2), 'Εξωτερικά Εκπεμπόμενη Ισχύς');
title('Εξωτερικά Εκπεμπόμενη Ισχύς');
%N
figure(15);
plotyy(t, Is, t, y(:,1)/Nth);
ax=plotyy(t, Is, t, y(:,1)/Nth);
xlabel('time (sec)');
ylabel(ax(1), 'I (current)');
ylabel(ax(2), 'Κανονικοποιημένος Αριθμός Φορέων');
title('Κανονικοποιημένος Αριθμός Φορέων');
%Eye
figure(16);
j=1;
hold on;
for i=1:1:31
   Pe_eye=[];
   while(t(j)<pulseDuration*i)
         Pe_eye=[Pe_eye, 0.5*vg*amir*h*f*y(j,2)];
         j=j+1;
   end
   plot(0:pulseDuration/length(Pe_eye):pulseDuration-pulseDuration/length(Pe_eye), Pe_eye);
end
xlabel('time (sec)');
ylabel('Εξωτερικά Εκπεμπόμενη Ισχύς');
title('Διάγραμμα Οφθαλμού Εκπεμπόνενης Ισχύος');

