%-------------------------------------
%%----- 1.b ---------------------------
%%
wavelength=1550*10^(-9);
n1eff_TE_1D=2.805103289;
n2eff_TE_1D=1.610716807;
n1eff_TM_1D=1.874758136;
n2eff_TM_1D=1.45;
c0=3*10^8;
m0=4*pi*10^(-7);
e0=8.854*10^(-12);
gamma_TE=zeros(1, 700+1);
gamma_TM=zeros(1, 700+1);
i=1;
for w=300*10^(-9):1*10^(-9):1000*10^(-9)
    x_all=linspace(-2*w, 3*w, 50000);
    x_slab=linspace(0, w, 10000);

    
    [ ~ , ~ , ~ , Ex_TM_all, Hy_TM_all ] = APDWG(wavelength , w , n1eff_TE_1D , n2eff_TE_1D , n2eff_TE_1D , x_all );
    P_slab=(1/2) * numerical_integration(Ex_TM_all(1,20001:30000).*(conj(Hy_TM_all(1,20001:30000))), x_all(20001:30000));
    P_all=(1/2) * numerical_integration(Ex_TM_all(1,:).*(conj(Hy_TM_all(1,:))), x_all);
    gamma_TE(i)=P_slab/P_all; 
    
    
    [ neff_TE_1 , ~ , Ey_TE_all , ~, ~ ] = APDWG(wavelength , w , n1eff_TM_1D , n2eff_TM_1D , n2eff_TM_1D , x_all );
    [ neff_TE_2 , ~ , Ey_TE_slab , ~, ~ ] = APDWG(wavelength , w , n1eff_TM_1D , n2eff_TM_1D , n2eff_TM_1D , x_slab );
    Hx_TM_all=-(neff_TE_1(1)/(c0*m0)).*Ey_TE_all(1,:);
    Hx_TM_slab=-(neff_TE_1(1)/(c0*m0)).*Ey_TE_slab(1,:);
    P_slab=(-1/2) * numerical_integration(Ey_TE_slab(1,:).*(conj(Hx_TM_slab)), x_slab);
    P_all= (-1/2) * numerical_integration(Ey_TE_all(1,:).*(conj(Hx_TM_all)), x_all);
    gamma_TM(i)=P_slab/P_all;
    
    i=i+1;
end
w=300*10^(-9):1*10^(-9):1000*10^(-9);
plot(w, gamma_TE);
hold on;
plot(w, gamma_TM);


%%------------------------------------
%%------- 1c--------------------------
%%
wavelength=1550*10^(-9);
nsi0=3.45;
sn_e=8.8*10^(-28);
sn_h=4.6*10^(-28);
nsi=@(N)(nsi0-sn_e*N-(sn_h*N)^(0.8));
w=500*10^(-9);
h=220*10^(-9);
t=50*10^(-9);
N_min=10^24;
N_max=10^26;
nsiO2=1.45;
n0=1;
i=1;
neff_TM0_2D=zeros(1, 100);
neff_TE0_2D=zeros(1, 100);

for N=(1:1:100)*10^24
    [ neff_TE_1D_1 , neff_TM_1D_1 , ~ , ~, ~ ] = APDWG( wavelength , h , nsi(N) , nsiO2 , n0 );
    [ neff_TE_1D_2 , neff_TM_1D_2 , ~ , ~, ~ ] = APDWG( wavelength , t , nsi(N) , nsiO2 , n0 );
    if(isempty(neff_TM_1D_2))
         neff_TM_1D_2=nsiO2;
    end
    [ ~ , neff_TM_2D , ~ , ~, ~ ] = APDWG( wavelength , w , neff_TE_1D_1(1) , neff_TE_1D_2(1) , neff_TE_1D_2(1) );
    neff_TM0_2D(i)=neff_TM_2D(1);
    [ neff_TE_2D , ~ , ~ , ~, ~ ] = APDWG( wavelength , w , neff_TM_1D_1(1) , neff_TM_1D_2(1) , neff_TM_1D_2(1) );
    neff_TE0_2D(i)=neff_TE_2D(1);
    
    i=i+1;
end
N=(1:1:100)*10^18;
semilogx(N, neff_TM0_2D);
hold on;
semilogx(N, neff_TE0_2D);
