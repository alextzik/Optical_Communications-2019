syms l v;
radius=4.1*10^(-6);
e0=8.854*10^(-12);
m0=4*pi*10^(-7);
a=[0.6961663 0.4079426 0.8974994];
b=[0.004629148 0.01351206 97.934062];
n=sqrt(1+a(1)*(l^2)/((l^2-b(1)))+a(2)*(l^2)/((l^2-b(2)))+a(3)*(l^2)/((l^2-b(3))));
c=1/sqrt(m0*e0);

%%
%--------------1.a-------------------
dn_2=diff(n,2);

Dm_array=zeros(1, 0.5/0.001);
i=1;
for lamdas=(1.2:0.001:1.7)*10^(-6);
    l=lamdas*10^6;
    Dm_array(i)=-lamdas*subs(dn_2)/c*10^18;
    i=i+1;
end
figure(1);
lamdas=(1.2:0.001:1.7)*10^(-6);
plot(lamdas, Dm_array);

%%
%-------------1.b---------------------
%Constants, Variables & Equations
delta=0.003;
k0=@(wl, index)(2*pi/(wl));
V=@(wl, index)(k0(wl, index)*radius*index*sqrt(2*delta));
char_eq=@ (b, V)((sqrt(1-b)*besselj(1,V*sqrt(1-b))/besselj(0, V*sqrt(1-b)))-(sqrt(b)*besselk(1,V*sqrt(b))/besselk(0, V*sqrt(b))));


%Find range of V for the given range of wavelengths
l=1.2;
Vmax=V(1200*10^(-9), double(subs(n)));
l=1.7;
Vmin=V(1700*10^(-9), double(subs(n)));
betas=zeros(1,int64(1.5/0.0001));

%Find numerically the data b=b(V)
i=1;
for Vs=1.5:0.0001:3
   characteristic_equation=@(b)(char_eq(b, Vs));
   betas(i)=abs(fsolve(characteristic_equation, 0.97)); 
   i=i+1; 
end
Vs=1.5:0.0001:3;

%{
%Curve fitting of the data
b_of_V=@(c,x)(c(1)*(c(4)*exp(x))./(c(5)+c(2)*exp(c(3)*x)));
c0=[1 0 1 0 1];
constants=nlinfit(Vs, betas, b_of_V, c0);
figure(2);
plot(Vs, betas);
hold on;
plot(Vs, b_of_V(constants,Vs));

%Find Dw
Vb=(v*b_of_V(constants, v));
dVb_2=diff(Vb,2);
%}

Dw_array=zeros(1, 0.5/0.001);
i=1;

%Other way
Vb_2=betas.*Vs;
dVb_1=gradient(Vb_2, 0.0001);
dVb_2=gradient(dVb_1, 0.0001);

for wl=(1.2:0.001:1.7)*10^(-6)
    l=wl*10^6;
    index=double(subs(n));
    v=V(wl, index);
    %Dw_array(i)=(-index*delta*subs(v)*double(subs(dVb_2))/(wl*c))*10^6;
    Dw_array(i)=(-index*delta*subs(v)*dVb_2(int64((roundn(v,-4)-1.5)*10000+1))/(wl*c))*10^6;
    i=i+1;
end
wl=(1.2:0.001:1.7)*10^(-6);
figure(2)
plot(wl, Dw_array);

%%
%-------------1.c-------------------
Dt_array=Dm_array+Dw_array;
figure(4)
plot(wl, Dt_array);

%%
%-------------1.4-------------------
%Get the function â(ù) from b=(â^2/k0-n2^2)/(n1^2-n2^2)  and b=b(V), V=V(ù)
l=1.7;
omega_min=Vmin/(sqrt(m0*e0)*radius*double(subs(n))*sqrt(2*delta));
l=1.2;
omega_max=Vmax/(sqrt(m0*e0)*radius*double(subs(n))*sqrt(2*delta));
propConsts=zeros(1,int64((omega_max*10^(-15)-omega_min*10^(-15)+0.02)/0.001));
i=1;

for omegas=(omega_min*10^(-15)-0.01:0.001:omega_max*10^(-15)+0.01)*10^15;
    wl=2*pi*c/omegas;
    l=wl*10^6;
    index=double(subs(n));
    v=V(wl,index);
    characteristic_equation=@(b)(char_eq(b, v));
    beta=abs(fsolve(characteristic_equation, 0.97)); 
    propConsts(i)=omegas/c*index*sqrt(1-2*delta*(1-beta));
    i=i+1;
end
omegas=(omega_min*10^(-15)-0.01:0.001:omega_max*10^(-15)+0.01)*10^15;

%%
%Get the necessary second derivative
dPropConsts_1=gradient(propConsts, 0.001);
dPropConsts_2=gradient(dPropConsts_1, 0.001);
DtFormula_array=zeros(1,0.5/0.001);
i=1;

for wl=(1.2:0.001:1.7)*10^(-6)
    omega=2*pi*c/wl;
    DtFormula_array(i)=-omega/wl*dPropConsts_2(roundn(omega*10^(-15)-(omega_min*10^(-15)-0.01),-3)*1000+1)*10^(-24);
    i=i+1;
end

wl=(1.2:0.001:1.7)*10^(-6);
figure(3);
hold on;
plot(wl, DtFormula_array); 
%%
%Curve Fit Dt
DtFormula_Array_Fit=(-1.628*10^(-6).*wl.^(-1.309)+82.23);
Dt_Array_fitted=(-2.573*10^(-6).*wl.^(-1.278)+84.55);

%------------------1.e-----------------
%%
S_1=gradient(Dt_Array_fitted, 0.001*10^(-6))*10^(-9);
S_2=gradient(DtFormula_Array_Fit, 0.001*10^(-6))*10^(-9);
figure(5);
wl=(1.2:0.001:1.7)*10^(-6);
plot(wl, S_1);
hold on;
wl=(1.2:0.001:1.7)*10^(-6);
plot(wl, S_2);