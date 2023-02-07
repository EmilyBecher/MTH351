%clear
%clc
clf

fontsize = 12;
%n = input('Give the number n defining the (2n + 1) degree Taylor Polynomial: ');
n = 13;
x = linspace(0,4*pi,1000);

%Calculate functions to be plotted
ym = sin(x);
[y0, y1] = eval_sin(n, x);

%Calculate the error of both approximations
err0 = max(abs(y0-ym));
err1 = max(abs(y1-ym));

%Generate the plot
figure(1)
axis = subplot(1,1,1);
set(axis, 'FontSize', 8);
ylim([-2, 2]);
hold on
plot(x,ym,'-K','LineWidth', 2, 'DisplayName', 'Matlab Sine')
plot(x,y0,'LineWidth', 2, 'DisplayName', strcat('Taylor Series, Max Error = ', num2str(err0, '%1.2e')))
plot(x,y1,'LineWidth', 2, 'DisplayName', strcat('Identities, Max Error = ', num2str(err1, '%1.2e')))
legend('Interpreter','latex','FontSize',fontsize,'Location','southwest')
title('Emily Becher', 'Plots of $y=sin(x)$ Using Various Approximations','Interpreter','latex','FontSize',fontsize)

%Calculates the Taylor coefficients
function [y0, y1] = eval_sin(n, x)
    coeff = ones(n+1, 1);
    sign = 1;
    fact = 1;

    coeff(1) = 1;
    for i = 2:n+1
        sign = -sign;
        fact = fact*(2*i-1)*(2*i-2);
        coeff(i) = sign/fact;
    end

    x1 = mod(x, 2*pi);
    x1 = x1-pi;
    y1 = -eval_poly(coeff, x1, n);
    y0 = eval_poly(coeff, x, n);
end

%Evaluates odd polynomials given coefficients
function value = eval_poly(coeff, x, n)

    xsq = x.*x;
    value = coeff(n+1)*ones(size(x));

    for i = n:-1:1
        value = coeff(i) + xsq.*value;
    end
    
    value = x.*value;

end