%clear
%clc
clf

fontsize = 12;
%n = input('Give the number n defining the degree of interpolation (uses n+1 input points): ');
n = 10;
x = linspace(-5,5,1000);
x1 = linspace(-5,5,n+1);
x2 = 5*cos(linspace(0,pi,n+1));

%Calculate functions to be plotted
ym = (1+x.^2).^(-1);
y0 = interpolate(n, x1);
y0 = flip(y0); %Reorder coefficients to use polyval
y0 = polyval(y0, x);
y1 = interpolate(n, x2);
y1 = flip(y1); %Reorder coefficients to use polyval
y1 = polyval(y1, x);

%Generate the plot
figure(1)
axis = subplot(1,1,1);
set(axis, 'FontSize', 8);
ylim([-1, 2]);
hold on
plot(x,ym,'-K','LineWidth', 2, 'DisplayName', 'Matlab Function')
plot(x,y0,'LineWidth', 2, 'DisplayName', 'Evenly Space Input Points')
plot(x,y1,'LineWidth', 2, 'DisplayName', 'Weighted Input Points')
legend({'Matlab Function','Evenly Space Input Points', 'Weighted Input Points'},'Location','northeast')
title('Emily Becher', 'Plots of $y=\frac{1}{1+x^2}$ Using Interpolation Approximations','Interpreter','latex','FontSize',fontsize)

%Calculates the coefficients
function y0 = interpolate(n, x)
    y = (1+x.^2).^(-1);
    y = transpose(y);
    A = zeros(n+1,n+1);
    for i = 1:1:n+1
        for j = 1:1:n+1
            A(i,j) = x(i)^(j-1);
        end
    end
    y0 = (A^-1)*y;
end

