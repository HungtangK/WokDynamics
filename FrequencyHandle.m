clear
close all
clc

load handle.mat

FHx_freq=FFT(FHx,t)
FHy_freq=FFT(FHy,t)
xH_freq=FFT(xH,t)
yH_freq=FFT(yH,t)

function Return=FFT(X,t)
    
    % Sampling frequency and length
    Fs=1/(t(2)-t(1));
    L=length(t);
    f = Fs*(0:(L/2))/L;

    % FFT
    Y = fft(X);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

    figure
    plot(f,P1)
    xlim([-0.1,10])

    % Max values
    A_max1=max(P1);
    f_max1=f(find(P1==A_max1));
    phi_max1=angle(Y(find(P1==A_max1)));
    P1(find(P1==A_max1))=0;


    A_max2=max(P1);
    f_max2=f(find(P1==A_max2));
    phi_max2=angle(Y(find(P1==A_max2)));
    P1(find(P1==A_max2))=0;

    A_max3=max(P1);
    f_max3=f(find(P1==A_max3));
    phi_max3=angle(Y(find(P1==A_max3)));
    P1(find(P1==A_max3))=0;

    % Plotting
    figure
    plot(t,X,'y','LineWidth',2)
    hold on 
    plot(t,A_max1*cos(2*pi*f_max1*t+phi_max1)+...
           A_max2*cos(2*pi*f_max2*t+phi_max2)+...
           A_max3*cos(2*pi*f_max3*t+phi_max3),'--r','LineWidth',2)
   legend('original','FFT')
    
    % return 
    Return = [ A_max1,f_max1,phi_max1,
            A_max2,f_max2,phi_max2,
            A_max3,f_max3,phi_max3];
end


