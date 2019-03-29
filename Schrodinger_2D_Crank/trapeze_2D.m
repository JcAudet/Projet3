function I=trapeze_2D(f,a,b,c,d,n,m)


if (n <= 0)
  error('n doit etre positif');
end

h = (b - a) / n;
k = (d - c) / m;
l = size(f);

if (l(2) ~= n+1 || l(1)~=m+1)
error('le nombre de points fournis ne correspond pas avec le nombre de sous-intervalles');
end

sum_abcd = f(1,1) + f(1,n+1) + f(m+1,1) + f(m+1,n+1);
sum_a = sum(f(2:m,1));
sum_b = sum(f(2:m,n+1));
sum_c = sum(f(1,2:n));
sum_d = sum(f(m+1,2:n));
F = f(2:m,2:n);

I = ( h * k / 4) * ( sum_abcd + 2*(sum_a + sum_b + sum_c + sum_d)+ 4*sum(sum(F)) );
