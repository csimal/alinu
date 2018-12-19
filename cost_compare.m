clc
clear
close all

f1 = @(m)( 8*m.^3 );
f2 = @(m)( 8*m.^3 - 4*m.^2 );
f3 = @(m)( 7*m.^3 );
f4 = @(m)( 7*m.^3 + 11*m.^2 );

m=1:30;

figure
subplot(2,1,1)
plot(m,f1(m))
hold on
plot(m,f3(m))
hold off
title('Product costs')

subplot(2,1,2)
plot(m,f2(m))
hold on
plot(m,f4(m))
hold off
title('Sum costs')