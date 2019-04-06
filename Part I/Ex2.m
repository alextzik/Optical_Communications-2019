%%-------------2.c----------------------
%%
char_eq=@ (b, V)((sqrt(1-b)*besselj(1,V*sqrt(1-b))/besselj(0, V*sqrt(1-b)))-(sqrt(b)*besselk(1,V*sqrt(b))/besselk(0, V*sqrt(b))));
betas=zeros(1,(12-0.1)/0.1);
i=1;
for V=0.1:0.1:12
   characteristic_equation=@(b)(char_eq(b, V));
   betas(i)=fsolve(characteristic_equation, 0.97); 
   i=i+1; 
end
V=0.1:0.1:12;
plot(V,abs(betas));

%%--------2.d---------------------------
%%
a=5*10^(-6);
w0=@(V)(a*(0.65+1.619*V^(-3/2)+2.879*V^(-6)));
eq=@(x, y, V)(exp((-(x^2+y^2)/(w0(V)^2))));
Aeff=zeros(1,1.6/0.02+1);
i=1;
for V=0.8:0.02:2.4
   disp(V);
   integrant_num=@(x, y)(abs(eq(x, y, V))^2);
   integrant_den=@(x, y)(abs(eq(x, y, V))^4);
   numerator=(integrate(integrant_num))^2;
   denominator=(integrate(integrant_den));
   Aeff(i)=numerator/denominator;
   i=i+1;
end
Aeff=Aeff./(pi*a^2);
V=0.8:0.02:2.4;
plot(V, Aeff);