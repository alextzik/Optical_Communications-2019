wavelength=1550*10^(-9);
Bt=10*10^9;
%%
%--------------------2.a-----------------------
%Graph the input gaussian pulse
b2=-2.1*10^(-26);
b3=0;
C=0;
A0=0.4;
T0=0.25*10^(-10);
FWHM=2*sqrt(log(2))*T0;
t=-4*FWHM:T0/100:4*FWHM;
pulseIn=A0*exp(-(1+1i*C)/2*(t/T0).^2);
figure(1);
plot(t, abs(pulseIn));


%FFT of input pulse
Fmax=1/(FWHM); %Considering that time duration and frequency spread are inversely related
Fs=10000*Fmax;
t=-4*FWHM:1/Fs:4*FWHM;
pulseIn=A0*exp(-(1+1i*C)/2*(t/T0).^2);
N=length(t);
FpulseIn=fft(pulseIn);
k=(0:1:N-1)-floor(N/2);
omegas=2*pi.*k.*Fs./N;
Max=max(abs(FpulseIn));
figure(2);

%Check if spectrum with fft is the same as analytic solution(set w=-w)
%omegas2=omegas(1:floor(length(omegas)/2));
analyticFpulseIn=A0*sqrt(2*pi*T0^2/(1+1i*C))*exp(-((omegas-2*pi).^2)*T0^2/(2*(1+1i*C)));
MaxA=max(abs(analyticFpulseIn));
hold on;
plot(omegas, abs(analyticFpulseIn));
energy=sum(abs(pulseIn).^2);
plot(omegas, fftshift(abs(FpulseIn)*sqrt(sum(abs(analyticFpulseIn).^2)/sum(abs(FpulseIn).^2))));

%%
%Check if ifft(FpulseIn)=pulseIn
%{
IFpulseIn=ifft(FpulseIn);

figure(1);
hold on;
plot(t, IFpulseIn);
%}

%Find spectrum of pulse in a given z using fft and analytic form
z=40000;
FpulseOut=fftshift(FpulseIn).*exp(1i*(0.5*b2*omegas.^2)*z);
pulseOut=(ifft(ifftshift(FpulseOut)));
figure(1);
hold on;
plot(t, abs(pulseOut));

t=-4*FWHM:1/Fs:4*FWHM;
analyticPulseOut=A0*T0/sqrt(T0^2-1i*b2*(1+1i*C)*z)*exp((-t.^2)*(1+1i*C)/(2*(T0^2-1i*b2*(1+1i*C)*z)));
plot(t, abs(analyticPulseOut)); 


%%
%---------------------2.b----------------------------
%Variables
pulseDuration=1/Bt;
b2=-2.1*10^(-26);
b3=1.3*10^(-40);
a=0.18/4.343;
nofPRBS=63;
prbs=[1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0];
Fs=120*10^9; %Since the spectrum [0,30]GHz is required

%Pulse Train Creation & Fourier Transform (2.b)
t = 0:1/Fs:nofPRBS*pulseDuration;
d = [0:pulseDuration:(nofPRBS-1)*pulseDuration; prbs]';
x = @(var)rectpuls(var-pulseDuration/2, pulseDuration);
pulseIn_b = pulstran(t,d,x);
figure(4);
plot(t, pulseIn_b);

figure(5);
FPulseIn_b=fft(pulseIn_b);
N=length(t);
k=(0:1:N-1)-floor(N/2);
omegas=2*pi.*k.*Fs./N;
energy= sum(pulseIn_b.^2);
plot(omegas(floor(N/2)+1:floor(N/2)+floor(N/4))./(2*pi), abs(FPulseIn_b(1:floor(N/4))));

%Optical Power
z=40000;
FpulseOut_2_40km=ifftshift(fftshift(FPulseIn_b).*exp(1i*(0.5*b2*(omegas.^2)-b3.*(omegas.^3)/6)*z));
pulseOut_2_40km=ifft(FpulseOut_2_40km);
power_2_40km=abs(pulseOut_2_40km).^2.*exp(-a*z/1000);
z=80000;
FpulseOut_2_80km=ifftshift(fftshift(FPulseIn_b).*exp(1i*(0.5*b2*omegas.^2-b3.*(omegas.^3)/6)*z));
pulseOut_2_80km=ifft(FpulseOut_2_80km);
power_2_80km=abs(pulseOut_2_80km).^2.*exp(-a*z/1000);
figure(6);
plot(t, power_2_40km);
hold on; 
plot(t, power_2_80km);
%plot(t, pulseIn_b);

%Eye diagrams
pulseOut_2_40km=sqrt(power_2_40km);
pulseOut_2_40km_eye = vec2mat(pulseOut_2_40km,pulseDuration*Fs);
figure(7);
plot(0:1/Fs:pulseDuration-1/Fs, pulseOut_2_40km_eye);

pulseOut_2_80km=sqrt(power_2_80km);
pulseOut_2_80km_eye = vec2mat(pulseOut_2_80km,pulseDuration*Fs);
figure(8);
plot(0:1/Fs:pulseDuration-1/Fs, pulseOut_2_80km_eye);


%%
%------------------------2.c----------------------------------------
%Variables
pulseDuration=1/Bt;
b2=-2.1*10^(-26);
b3=1.3*10^(-40);
a=0.18/4.343;
nofPRBS=63;
prbs=[1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0];
Fs=120*10^9; %Since the spectrum [0,30]GHz is required

%Pulse Train Creation & Fourier Transform (2.c)
t = 0:1/Fs:nofPRBS*pulseDuration;
d = [0:pulseDuration:(nofPRBS-1)*pulseDuration; prbs]';
x = @(var)gausPulse(var, pulseDuration, 0, 20*10^(-12), 1);
pulseIn_c = pulstran(t,d,x);
figure(9);
plot(t, pulseIn_c);
figure(10);
FPulseIn_c=fft(pulseIn_c);
N=length(t);
k=(0:1:N-1)-floor(N/2);
omegas=2*pi.*k.*Fs./N;
energy= sum(pulseIn_c.^2);
plot(omegas(floor(N/2)+1:floor(N/2)+floor(N/4))./(2*pi), abs(FPulseIn_c(1:floor(N/4))));

%Optical Power
z=40000;
FpulseOut_3_40km=ifftshift(fftshift(FPulseIn_c).*exp(1i*(0.5*b2*omegas.^2-b3.*(omegas.^3)/6)*z));
pulseOut_3_40km=ifft(FpulseOut_3_40km);
power_3_40km=abs(pulseOut_3_40km).^2.*exp(-a*z/1000);
z=80000;
FpulseOut_3_80km=ifftshift(fftshift(FPulseIn_c).*exp(1i*(0.5*b2*omegas.^2-b3.*(omegas.^3)/6)*z));
pulseOut_3_80km=ifft(FpulseOut_3_80km);
power_3_80km=abs(pulseOut_3_80km).^2.*exp(-a*z/1000);
figure(11);
plot(t, power_3_40km);
hold on; 
plot(t, power_3_80km);


%Eye diagrams
pulseOut_3_40km=sqrt(power_3_40km);
pulseOut_3_40km_eye = vec2mat(pulseOut_3_40km,pulseDuration*Fs);
figure(12);
plot(0:1/Fs:pulseDuration-1/Fs, pulseOut_3_40km_eye);

pulseOut_3_80km=sqrt(power_3_80km);
pulseOut_3_80km_eye = vec2mat(pulseOut_3_80km,pulseDuration*Fs);
figure(13);
plot(0:1/Fs:pulseDuration-1/Fs, pulseOut_3_80km_eye);


%%
%---------------------2.d-------------------------------------
%Variables
pulseDuration=1/Bt;
b2=-2.1*10^(-26);
b3=1.3*10^(-40);
a=0.18/4.343;
nofPRBS=63;
prbs=[1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0];
Fs=120*10^9; %Since the spectrum [0,30]GHz is required

%Find optimal To
x=0:0.001*10^(-12):100*10^(-12);
To40s=abs(x).*sqrt(1+(40000./((x.^2)./abs(b2))).^2);
[m, j]=min(To40s);
bestTo_40km=x(j);

To80s=((abs(x).*sqrt(1+(80000./((x.^2)./abs(b2))).^2)));
[m, j]=min(To80s);
bestTo_80km=x(j);

%Calculations with optimal To_40km
%Pulse Train Creation & Fourier Transform (2.c)
t = 0:1/Fs:nofPRBS*pulseDuration;
d = [0:pulseDuration:(nofPRBS-1)*pulseDuration; prbs]';
x = @(var)gausPulse(var, pulseDuration, 0, bestTo_40km, 1);
pulseIn_d_40 = pulstran(t,d,x);
figure(14);
plot(t, pulseIn_d_40);
figure(15);
FPulseIn_d_40=fft(pulseIn_d_40);
N=length(t);
k=(0:1:N-1)-floor(N/2);
omegas=2*pi.*k.*Fs./N;
energy= sum(pulseIn_d_40.^2);
plot(omegas(floor(N/2)+1:floor(N/2)+floor(N/4))./(2*pi), abs(FPulseIn_d_40(1:floor(N/4))));

%Optical Power
z=40000;
FpulseOut_4_1_40km=ifftshift(fftshift(FPulseIn_d_40).*exp(1i*(0.5*b2*omegas.^2-b3.*(omegas.^3)/6)*z));
pulseOut_4_1_40km=ifft(FpulseOut_4_1_40km);
power_4_1_40km=abs(pulseOut_4_1_40km).^2.*exp(-a*z/1000);
z=80000;
FpulseOut_4_1_80km=ifftshift(fftshift(FPulseIn_d_40).*exp(1i*(0.5*b2*omegas.^2-b3.*(omegas.^3)/6)*z));
pulseOut_4_1_80km=ifft(FpulseOut_4_1_80km);
power_4_1_80km=abs(pulseOut_4_1_80km).^2.*exp(-a*z/1000);
figure(16);
plot(t, power_4_1_40km);
hold on; 
plot(t, power_4_1_80km);

%Eye diagrams
pulseOut_4_1_40km=sqrt(power_4_1_40km);
pulseOut_4_1_40km_eye = vec2mat(pulseOut_4_1_40km,pulseDuration*Fs);
figure(17);
plot(0:1/Fs:pulseDuration-1/Fs, pulseOut_4_1_40km_eye);

pulseOut_4_1_80km=sqrt(power_4_1_80km);
pulseOut_4_1_80km_eye = vec2mat(pulseOut_4_1_80km,pulseDuration*Fs);
figure(18);
plot(0:1/Fs:pulseDuration-1/Fs, pulseOut_4_1_80km_eye);

%%
%Calculations with optimal To_80km
%Pulse Train Creation & Fourier Transform (2.c)
t = 0:1/Fs:nofPRBS*pulseDuration;
d = [0:pulseDuration:(nofPRBS-1)*pulseDuration; prbs]';
x = @(var)gausPulse(var, pulseDuration, 0, bestTo_80km, 1);
pulseIn_d_80 = pulstran(t,d,x);
figure(19);
plot(t, pulseIn_d_80);
figure(20);
FPulseIn_d_80=fft(pulseIn_d_80);
N=length(t);
k=(0:1:N-1)-floor(N/2);
omegas=2*pi.*k.*Fs./N;
energy= sum(pulseIn_d_80.^2);
plot(omegas(floor(N/2)+1:floor(N/2)+floor(N/4))./(2*pi), abs(FPulseIn_d_80(1:floor(N/4))));

%Optical Power
z=40000;
FpulseOut_4_2_40km=ifftshift(fftshift(FPulseIn_d_80).*exp(1i*(0.5*b2*omegas.^2-b3.*(omegas.^3)/6)*z));
pulseOut_4_2_40km=ifft(FpulseOut_4_2_40km);
power_4_2_40km=abs(pulseOut_4_2_40km).^2.*exp(-a*z/1000);
z=80000;
FpulseOut_4_2_80km=ifftshift(fftshift(FPulseIn_d_80).*exp(1i*(0.5*b2*omegas.^2-b3.*(omegas.^3)/6)*z));
pulseOut_4_2_80km=ifft(FpulseOut_4_2_80km);
power_4_2_80km=abs(pulseOut_4_2_80km).^2.*exp(-a*z/1000);
figure(21);
plot(t, power_4_2_40km);
hold on; 
plot(t, power_4_2_80km);


%Eye diagrams
pulseOut_4_2_40km=sqrt(power_4_2_40km);
pulseOut_4_2_40km_eye = vec2mat(pulseOut_4_2_40km,pulseDuration*Fs);
figure(22);
plot(0:1/Fs:pulseDuration-1/Fs, pulseOut_4_2_40km_eye);

pulseOut_4_2_80km=sqrt(power_4_2_80km);
pulseOut_4_2_80km_eye = vec2mat(pulseOut_4_2_80km,pulseDuration*Fs);
figure(23);
plot(0:1/Fs:pulseDuration-1/Fs, pulseOut_4_2_80km_eye);


%%
%-------------------------2.e-----------------------------------------
%Variables
pulseDuration=1/Bt;
b2=-2.1*10^(-26);
b3=1.3*10^(-37);
a=0.18/4.343;
nofPRBS=63;
prbs=[1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0];
Fs=120*10^9; %Since the spectrum [0,30]GHz is required

%Pulse Train Creation & Fourier Transform (2.c)
t = 0:1/Fs:nofPRBS*pulseDuration;
d = [0:pulseDuration:(nofPRBS-1)*pulseDuration; prbs]';
x = @(var)gausPulse(var, pulseDuration, 0, 20*10^(-12), 1);
pulseIn_c = pulstran(t,d,x);
FPulseIn_c=fft(pulseIn_c);

%Optical Power
z=80000;
FpulseOut_3_80km=ifftshift(fftshift(FPulseIn_c).*exp(1i*(0.5*b2*omegas.^2-b3.*(omegas.^3)/6)*z));
pulseOut_3_80km=ifft(FpulseOut_3_80km);
power_3_80km=abs(pulseOut_3_80km).^2.*exp(-a*z/1000);

%Eye diagram
pulseOut_3_80km=sqrt(power_3_80km);
pulseOut_3_80km_eye = vec2mat(pulseOut_3_80km,pulseDuration*Fs);
figure(37);
plot(0:1/Fs:pulseDuration-1/Fs, pulseOut_3_80km_eye);

%%Biggest change from b3=1.3*10^(-37)->b3=1.3*10^(-36)