clear
clc
clf

fontsize = 12;
%n = input('Give the number n defining the degree of interpolation (uses n+1 input points): ');
m = 400;
n = 10^9;
a = 0;
b = 1;
x = linspace(0,1,1000);
exact = (exp(1)-1)/m;
num = 1:1:log(n);
for i = 1:1:log(n)
    num(i) = 2^i;
end

z = length(num);
%Presize vectors
err_mid = zeros(1,z);
err_trap = zeros(1,z);
err_simp = zeros(1,z);

for j = 1:1:z
    h = (b-a)/(2^j);
    %Compute Nodes
    x_k = linspace(a,b,2^j+1);
    t_k = zeros(1,2^j);
    for i = 1:1:2^j
        t_k(i) = x_k(i) + h/2;
    end
    
    %Compute y vectors
    y_t = t_k.^(m-1).*exp(t_k.^m);
    y_x = x_k.^(m-1).*exp(x_k.^m);
    
    %Calculate Integrals
    I_mid = midpoint(y_t, h, 2^j);
    I_trap = trapezoid(y_x, h, 2^j);
    I_simp = simpson(y_t, y_x, h, 2^j);

    %Calculate error
    err_mid(j) = abs(exact-I_mid);
    err_trap(j) = abs(exact-I_trap);
    err_simp(j) = abs(exact-I_simp);
end

%Generate the plot
figure(1)
axis = subplot(1,1,1);
set(axis, 'FontSize', 8);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
%ylim([0, 1]);
%xlim([-1,n]);
hold on
plot(num, err_mid, 'o','LineWidth', 2)
plot(num, err_trap, 'o','LineWidth', 2)
plot(num, err_simp, 'o','LineWidth', 2)
legend({'Midpoint', 'Trapezoid', 'Simpson'},'Location','northeast')
xlabel('Number of Function Approximations','Interpreter','latex')
ylabel('Error','Interpreter','latex')
title('Emily Becher', 'Error in integrating $f(x) = x^{400-1}e^{x^{400}}$','Interpreter','latex','FontSize',fontsize)
%print -deps exp_convergence

%Approximates integrals using midpoint rule
function s_mid = midpoint(y_t, h, n)
    s_mid = 0;
    for i = 1:1:n
        s_mid = s_mid + y_t(i);
    end
    s_mid = h*s_mid;
end

%Approximates integrals using trapezoid rule
function s_trap = trapezoid(y_x, h, n)
    s_trap = y_x(1) +y_x(n+1);
    for i = 2:1:n
        s_trap = s_trap + 2*y_x(i);
    end
    s_trap = h*s_trap/2;
end

%Approximates integrals using Simpson's rule
function s_simp = simpson(y_t, y_x, h, n)
   s_simp = y_x(1) + y_x(n+1); 
   for i = 1:1:n
       s_simp = s_simp + 4*y_t(i);
   end
   for i = 2:1:n
       s_simp = s_simp + 2*y_x(i);
   end
   s_simp = h*s_simp/6;
end



