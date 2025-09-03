function Task1_Bobrovnikova
  clear;clc;

  C = -1/4; D = -1; N = 10;
  xmin = -C; xmax = 100;
  ymin = -0.5; ymax = 0;
  deltax = deltay = 1e-2;

  t = xmin:deltax:xmax;
  plot(t, D./(t+C), '--k', 'DisplayName' , 'The correct curve');
  hold on

  t = randi([ceil(-C) 100],N,1);
  x = sort(t);
  y = D./(x+C);
  Eps = 0.2*abs(y(N));
  noise = -Eps + 2.*Eps.*rand(N,1);
  y += noise;

  A = [y -y.^0];
  [Q R] = qr(A, 0);
  u = R\Q'*(-x.*y);
  C = u(1)
  D = u(2)

  [X,Y] = meshgrid(xmin:deltax:xmax,ymin:deltay:ymax);
  Z = X.*Y + C*Y - D;
  contour(X,Y,Z,[0 0], 'b', 'LineWidth', 2, 'DisplayName' , 'The resulting curve')
  c.LineWidth = 2;
  hold on
  plot(x,y,'pentagramk', 'DisplayName' , 'Starting points with noise');
  title( 'xy + Cy = D' )
  xlabel('X');
  ylabel('Y');
  legend('Location', 'southeast', 'FontSize' , 12)
  legend( 'boxoff' )
  end
