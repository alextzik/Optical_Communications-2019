function I = numerical_integration(y,x)
I = 0;
dx = x(2)-x(1);
for i = 1:1:length(x)
    I = I + y(i)*dx;
end

end