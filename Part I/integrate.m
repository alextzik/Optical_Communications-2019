function [I] = integrate(eq)
    I=0;
    dx=0.5*10^(-6);
    dy=0.5*10^(-6);
    for x=(-100:0.5:100)*10^(-6)
       for y=(-100:0.5:100)*10^(-6) 
            I=I+eq(x,y)*dx*dy;
       end 
    end
end