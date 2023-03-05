%clear
%clc
clf

fontsize = 12;
%n = input('Give the number n defining the degree of interpolation (uses n+1 input points): ');
n = 5;
nodes = 1:1:n;
x = linspace(-1,1,500);
exact = abs(x);
err = ones(1,n);
err_c = ones(1,n);

for k = 1:n
    %Compute Nodes
    X = linspace(-1,1,k) 
    ints = 0:1:k-1;
    X_c = cos((pi/2)*(2*ints+1)/(k)); %Chebyshev spacing
    Y = abs(X);
    Y_c = abs(X_c);
    
    %Calculate divided differences
    D = div_diff(X, Y);
    D_c = div_diff(X_c, Y_c);
    
    %Inerpolate divided differences
    y = interp(x, X, D);
    y_c = interp(x, X_c, D_c);
    
    %Calculate error
    err(k) = max(abs(exact-y));
    err_c(k) = max(abs(exact-y_c));
end

%Generate the plot
figure(1)
axis = subplot(1,1,1);
set(axis, 'FontSize', 8);
%set(gca, 'YScale', 'log');
%ylim([0, 1]);
%xlim([-1,n]);
hold on
%plot(nodes,err, 'o','LineWidth', 2)
%plot(nodes,err_c, 'o', 'LineWidth', 2)
plot(x, exact, '-k');
plot(x, y);
plot(x, y_c);
legend({'Even Spacing', 'Chebyshev'},'Location','northeast')
xlabel('Number of Nodes $(n+1)$','Interpreter','latex')
ylabel('$||f-P_n||_\infty$','Interpreter','latex')
title('Emily Becher', 'Error in interpolating polynomials for $f(x) = e^x$','Interpreter','latex','FontSize',fontsize)
%print -deps exp_convergence

%Calculates the divided differences
function d = div_diff(x, y)
    n = length(x);
    d = y;
    for i = 2:n
        for j = n:-1:i
            d(j) = (d(j)-d(j-1))/(x(j)-x(j-i+1));
        end
    end
end

%Newton's Divided Differences Interpolation
function y = interp(x, X, Y)
    n = length(X);
    y = Y(n)*ones(size(x));
    for i = n-1:-1:1
        y = Y(i) + (x - X(i)).*y;
    end
end

