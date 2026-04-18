[x,y,z] = ABM_complex_lorenz(10,28,8/3,0.995);
plot(real(x), real(y))
grid on
xlabel('Re(x)')
ylabel('Re(y)')
zlabel('z')