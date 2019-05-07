function [ output ] = gausPulse( t1, pulseDuration, C, T0, A0 )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    t=t1+pulseDuration/2;
    output=zeros(1, length(t));
    for i=1:1:length(t)
        if t(i)<pulseDuration
            output(i)=A0*exp(-(1+1i*C)/2*(t(i)/T0).^2);
        end
    end
end

