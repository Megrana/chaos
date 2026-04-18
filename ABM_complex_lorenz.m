function [x,y,z]=ABM_complex_lorenz(a1,a2,a3,q)

if nargin==0
    a1=10; a2=28; a3=8/3;
    q=0.995;
end

h=0.01;
N=5000;

% ---------- 初值（复数） ----------
x0 = 7 + 8i;
y0 = 5 + 6i;
z0 = 12;

% ---------- 预分配 ----------
x = zeros(1,N+1); 
y = zeros(1,N+1); 
z = zeros(1,N+1);

x1 = zeros(1,N+1); 
y1 = zeros(1,N+1); 
z1 = zeros(1,N+1);

% ---------- 定义系统 ----------
f1 = @(x,y,z) a1*(y - x);
f2 = @(x,y,z) a2*x - y - x*z;
f3 = @(x,y,z) 0.5*(x*conj(y) + conj(x)*y) - a3*z;

% ---------- 第一步（启动 ABM） ----------
x1(1)=x0 + h^q * f1(x0,y0,z0)/(gamma(q)*q);
y1(1)=y0 + h^q * f2(x0,y0,z0)/(gamma(q)*q);
z1(1)=z0 + h^q * f3(x0,y0,z0)/(gamma(q)*q);

x(1)=x0 + h^q*(f1(x1(1),y1(1),z1(1)) + q*f1(x0,y0,z0))/gamma(q+2);
y(1)=y0 + h^q*(f2(x1(1),y1(1),z1(1)) + q*f2(x0,y0,z0))/gamma(q+2);
z(1)=z0 + h^q*(f3(x1(1),y1(1),z1(1)) + q*f3(x0,y0,z0))/gamma(q+2);

% ---------- 主循环 ----------
for n=1:N
    
    M1=0; M2=0; M3=0;
    N1=0; N2=0; N3=0;

    for j=1:n
        b = (n-j+1)^q - (n-j)^q;
        a = (n-j+2)^(q+1) + (n-j)^(q+1) - 2*(n-j+1)^(q+1);

        fx = f1(x(j),y(j),z(j));
        fy = f2(x(j),y(j),z(j));
        fz = f3(x(j),y(j),z(j));

        N1 = N1 + b*fx;
        N2 = N2 + b*fy;
        N3 = N3 + b*fz;

        M1 = M1 + a*fx;
        M2 = M2 + a*fy;
        M3 = M3 + a*fz;
    end

    % ---------- 预测 ----------
    x1(n+1)=x0 + h^q*N1/(gamma(q)*q);
    y1(n+1)=y0 + h^q*N2/(gamma(q)*q);
    z1(n+1)=z0 + h^q*N3/(gamma(q)*q);

    % ---------- 校正 ----------
    x(n+1)=x0 + h^q*(f1(x1(n+1),y1(n+1),z1(n+1)) + M1)/gamma(q+2);
    y(n+1)=y0 + h^q*(f2(x1(n+1),y1(n+1),z1(n+1)) + M2)/gamma(q+2);
    z(n+1)=z0 + h^q*(f3(x1(n+1),y1(n+1),z1(n+1)) + M3)/gamma(q+2);

end