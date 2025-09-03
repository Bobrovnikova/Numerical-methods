function Task2_Bobrovnikova
  clear;
  clc;
  %y(1) = x, y(2) = y, y(3) = v_x, y(4) = v_y
  %f(1) = v_x, f(2) = v_y, f(3) = a_x, f(4) = a_y
  r= @(y) sqrt(y(1)^2+y(2)^2); %|r|
  f = @(y)[y(3);y(4); -4*pi^2*y(1)/r(y)^3; -4*pi^2*y(2)/r(y)^3];
  L = @(y) y(1)*y(4) - y(2)*y(3); %векторное произведение r*v = 2*dS/dt,
  %равное двум секторным скоростям
  T = 0.5;
  h = 0.0001; %разбиение
  N = floor(T/h);
  %вектор точек по оси времени
  t = 0: h : T;
  y = zeros(4, N+1);
  z = zeros(1,N+1);
  u0 = [1 0 0 pi; 1 0 0 2.5; 1 0 0 4; 1.7 0 0 1.7; 3 0 0 1.8];%5 начальных условий
  y(:, 1) = u0(2, :)';
##  %метод Верле начало:
##  y(:,2)=Euler_step(f,y(:,1),h);%старт методом Эйлера
##  for (n = 2: N)
##    y(1, n + 1) = Verlet_step(y(1, n),y(2, n), y(1, n-1), h);
##    y(2, n + 1) = Verlet_step(y(2, n),y(1, n), y(2, n-1), h);
##  endfor
##  for (n = 2: N)%восстановление скоростей через разностные производные
##    y(3, n) = (y(1,n+1)-y(1,n-1))/(2*h);
##    y(4, n) = (y(2,n+1)-y(2,n-1))/(2*h);
##  endfor
##  y(3, N+1) = (y(1,N+1)-y(1,N))/(h);
##  y(4, N+1) = (y(2,N+1)-y(2,N))/(h);
##  %конец метода Верле
  %начало метода RK2:
  for (n = 1: N)
    y(:, n + 1) = RK2_step(f, y(:, n), h);
  endfor
  %конец метода RK2
  for (n = 1: N+1)
    z(n) = L(y(:, n));
endfor
y(1,N+1)
y(2,N+1)
  figure(1)%орбита
  hold on;
  plot (y(1, :), y(2, :), 'g',
        [0], [0], 'ro')
  grid on
  grid minor
  axis equal;
  title('Метод RK2, u0 = [3; 0; 0; 1.8]')
  xlabel('x')
  ylabel('y')
  hold off;

  figure(2)%секторная скорость
  hold on;
  plot (t, z, '.m')
  title('Метод Верле, u0 = [1; 0; 0; pi]')
  xlabel('t')
  ylabel('L = 2*dS/dt')
  %axis equal;
  hold off;
end;

function y_next = RK2_step(f, y, h)
  y_next = y + h * f(y + h*f(y)/2);
end

function x2=Verlet_step(x1,y1, x0, h)
  x2=2*x1-x0+h^2*-4*pi^2*x1/(sqrt(x1^2+y1^2))^3;
end

function y_next=Euler_step(f,y,h)
  y_next=y+h*f(y);
  end
