% close all;
clear;

N = 33000;
tt = linspace(0,21,N);
Fs = 1/(tt(2)-tt(1));

yy1 = 3.5*sin(123.9/30*2*pi*tt);
% yy1 = yy1.*(yy1>0);

yy2 = 2.0*sin(127.2/30*2*pi*tt);
% yy2 = yy2.*(yy2>0);

yy3 = 0.98*sin(120.3/30*2*pi*tt);
% yy3 = yy3.*(yy3>0);

yy4 = 0.4177*sin(3.99*2*pi*tt);


figure(5)
plot(tt,yy1+yy2+yy3)
title('Toy beatwave')

xlim([7 21])
ylim([0 inf])

ff = fft(yy1+yy2+yy3);

P2 = abs(ff/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
ff_th = ff(1:N/2+1);
ff_th(ff_th<max(ff_th)/100) = 0;
theta = angle(ff_th);
f = Fs*(0:(N/2))/N;
figure(6)
plot(f,P1)
title('FFT of toy model')
xlim([0 20])

[peak,loc] = max(P1.*(f>2));
freq = f(loc);
th = theta(loc);


figure(7)
hold on
plot(tt,peak*cos(2*pi*freq*tt + th))
plot(tt,yy1+yy2)
xlim([0 20])

hold off

figure(8)
plot(f,theta)
title('tetilla')
xlim([0 20])


if false
L = 1000;
t = linspace(0, 1, L);
f1 = 10;
f2 = 20;
f3 = 15;
s1 = sin(2*pi*f1*t);
s2 = cos(2*pi*f2*t);
s3 = sin(2*pi*f3*t) + cos(2*pi*f3*t);
Ts = mean(diff(t));
Fs = 1/Ts;
Fn = Fs/2;
Fs1 = fft(s1)/L;
Fs2 = fft(s2)/L;
Fs3 = fft(s3)/L;
Fv = linspace(0, 1, fix(L/2)+1)*Fn;
Iv = 1:length(Fv);
Fs12 = fft(s1+s2)/L;
figure(1)
subplot(2,1,1)
plot(Fv, abs(Fs12(Iv)), '-b')
axis([0  50    ylim])
grid
title('sin(f_1) + cos(f_2)')
subplot(2,1,2)
plot(Fv, angle(Fs12(Iv))*180/pi, '-b')
axis([0  50    -180  180])
grid
figure(2)
subplot(2,1,1)
plot(Fv, abs(Fs3(Iv)), '-g')
axis([0  50    ylim])
grid
title('sin(f_3) + cos(f_3)')
subplot(2,1,2)
plot(Fv, angle(Fs3(Iv))*180/pi, '-g')
axis([0  50    -180  180])
grid
figure(3)
subplot(4,1,1)
plot(Fv, abs(Fs1(Iv)), '-b')
axis([0  50    ylim])
grid
title('sin(f_1)')
subplot(4,1,2)
plot(Fv, angle(Fs1(Iv))*180/pi, '-b')
axis([0  50    -180  180])
grid
subplot(4,1,3)
plot(Fv, abs(Fs2(Iv)), '-r')
axis([0  50    ylim])
grid
title('cos(f_2)')
subplot(4,1,4)
plot(Fv, angle(Fs2(Iv))*180/pi, '-r')
axis([0  50    -180  180])
grid
end