clear
clc
clf

fontsize = 12;
%n = input('Give the number n defining the degree of interpolation (uses n+1 input points): ');
m = 500;
n = 15;
a = 0;
b = 1;
period = 1;
x = linspace(0,1,1000);
num = linspace(1,n,n);

%Square wave
exact_square = zeros(size(x));
for i = 1:1:500
    exact_square(i) = 1;
end

%polynomial
exact_poly = x.^2;

%Triangular wave
exact_tri = zeros(size(x));
for i = 1:1:500
    exact_tri(i) = 2*x(i);
end
for i = 500:1:1000
    exact_tri(i) = -2*x(i)+2;
end

%Calculate Fourier Series approximations
Fourier_square = FS_approx(x, exact_square, n, a, b, period);
Fourier_poly = FS_approx(x, exact_poly, n, a, b, period);
Fourier_tri = FS_approx(x, exact_tri, n, a, b, period);


%Generate the plot
figure(1)
axis = subplot(1,1,1);
set(axis, 'FontSize', 8);
hold on
plot(x, exact_square, '-k','LineWidth', 2)
plot(x, Fourier_square);
legend({'Square Wave', 'Fourier Series: n = 100'},'Location','northeast')
%xlabel('Number of Fourier Series Coefficients','Interpreter','latex')
%ylabel('Error','Interpreter','latex')
title('Fourier Series Approximations of a Square Wave','Interpreter','latex','FontSize',fontsize)

figure(2)
axis = subplot(1,1,1);
set(axis, 'FontSize', 8);
ylim([0,1.2]);
hold on
plot(x, exact_poly, '-k','LineWidth', 2)
plot(x, Fourier_poly);
legend({'$y=x^2$', 'Fourier Series: n = 15'},'Location','northeast','Interpreter','latex')
%xlabel('Number of Fourier Series Coefficients','Interpreter','latex')
%ylabel('Error','Interpreter','latex')
title('Fourier Series Approximation of $y=x^2$','Interpreter','latex','FontSize',fontsize)

figure(3)
axis = subplot(1,1,1);
set(axis, 'FontSize', 8);
hold on
plot(x, exact_tri, '-k','LineWidth', 2)
plot(x, Fourier_tri);
legend({'Triangular Wave', 'Fourier Series: n = 100'},'Location','northeast','Interpreter','latex')
%xlabel('Number of Fourier Series Coefficients','Interpreter','latex')
%ylabel('Error','Interpreter','latex')
title('Fourier Series Approximation of a Triangular Wave','Interpreter','latex','FontSize',fontsize)

%Approximates integrals using trapezoid rule
function s_trap = trapezoid(y_x, h)
    n = length(y_x) - 1;
    s_trap = y_x(1) +y_x(n+1);
    for i = 2:1:n
        s_trap = s_trap + 2*y_x(i);
    end
    s_trap = h*s_trap/2;
end

%Calculate average value
function Fourier = FS_approx(x, exact, n, a, b, T)  
    freq = 1/T;
    w = 2*pi()*freq;
    h = (b-a)/length(x);
    avg = trapezoid(exact, h);
    avg = avg/T;
    
    a_n = zeros(1,n);
    b_n = zeros(1,n);
    for i = 1:1:n
        %Calculate cos coeffs
        y_x = exact.*cos(w*i*x);
        a_n(i) = 2/T*trapezoid(y_x, h);
        
        %Calculate sin coeffs
        y_x = exact.*sin(w*i*x);
        b_n(i) = 2/T*trapezoid(y_x, h);
    end
    
    %Evaluate sinusoid
    Fourier = zeros(size(x));
    for i = 1:1:n
        Fourier = Fourier + a_n(i)*cos(w*i*x) + b_n(i)*sin(w*i*x);
    end
    Fourier = Fourier + avg;
end



