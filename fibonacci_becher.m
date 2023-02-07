function [ans1, ans2] = fibonacci_becher(n, d)
%First method
m = 10^d;
prev = 0;
fib = 1;
if (n == 0)
    fib = 0;
end
if (n == 1)
    fib = 1;
else
    for i = 1:n-1
        w = uint64(fib + prev);
        prev = fib;
        fib = mod(w, 10^17);
    end 
end
ans1 = mod(fib, m);

%Second Method
A = [0 1; 1 1];
B = [0; 1];
C = (A^(n-1))*B; %doubles
fib = C(2,1);
ans2 = mod(fib, m);