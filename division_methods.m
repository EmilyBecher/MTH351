%clear
%clc
clf

fontsize = 18;
x = input('Give the number z for the division 1/z: ');
%x = 400;
%n = input('Give the maximum number of iterations: ');
n = 50;
b_error = zeros(1,n+1);
n_error = zeros(1,n+1);

%Find interval [a,b]
z = double(x);
y = typecast(z, 'uint64');
w = bitget(y,63:-1:53);
v = double(bit2int(w',11));
r = double(v - 1023);
a = 0.5^(r+1);
b = 0.5^r;

%Calculate 1/z and absolute error
[b_root, b_error] = bisect(a, b, z, n, r);
display(b_root)
[n_root, n_error] = newton(a, z, n, r);
display(n_root)
e_root = 1/z

%Generate Plot
%x = linspace(1,n,n);
x = 1:1:n+1;
figure(1)
axis = subplot(1,1,1);
set(axis, 'FontSize', 8);
set(gca, 'YScale', 'log');
xlim([0, n+2]);
hold on
plot(x, b_error, '-o', 'LineWidth', 2, 'DisplayName', 'Bisection Method')
plot(x, n_error, '-o', 'LineWidth', 2, 'DisplayName', 'Newton Method')
legend({'Bisection Method','Newton-Raphson Method'},'Location','northeast')
xlabel('Number of Iterations')
ylabel('Absolute Error')
title('Emily Becher', strcat('Error in computing $1/z$ where $z=$ ', num2str(z)),'Interpreter','latex','FontSize',fontsize)

function [root1, error] = bisect(a, b, z, n, e)
    error = zeros(1,n+1);
    for i = 0:1:n
        fb = 1 - b*z;
        c = 0.5*(a+b);
        fc = 1 - z*c;
        if fb == 0
            root1 = b;
        else 
            if sign(fc)*sign(fb) < 0
                a = c;
            else
                b = c;
            end
            root1 = c;
        end
        error(i+1) = abs(1/z-root1);
        if error(i+1) == 0
            error(i+1) = 2^(-53-e);
        end
    end
    root1 = c;
end

function [root2, error] = newton(g, z, n, e)
    error = ones(1,n+1);
    for i = 0:1:n
        x = 2*g-z*g*g;
        g = x;
        error(i+1) = abs(1/z-g);
        if error(i+1) == 0
            error(i+1) = 2^(-53-e);
        end
    end    
    root2 = g;
end