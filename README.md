# ProyectoAN2
Rafael Calderón; Estely Méndez



% Nombres: -Rafael Aníbal Calderón Melara  CM11099
%          -Estely Dolores Méndez Moreno   MM09204

% Fecha de creación: 15 de Noviembre de 2016
% Esta función aproxima la solución de la ecuación de Poisson
%      f_xx(x,y) + f_yy(x,y) = f(x,y); con x en [a,b], y en [c,d]
% sujeto a las condiciones iniciales
%      u(x,y) = g(x,y) si x=a o x=b & y pertenece a [c,d]
%      u(x,y) = g(x,y) si y=c o y=d & x pertenece a [a,b]
% El algoritmo utiliza el método de diferencias finitas 
% y a su vez el método de Gauss Seidel.
% Notar que se utiliza la norma infinita en nuestro código, 
% pero se podría utilizar cualquier otra norma.
% INPUT:  -- f,g: Funciones
%        -- a,b,c,d: Valores frontera.
%        -- n>=3, m>=3: Enteros que sirven para determinar las 
%                   dimensiones de la cuadrícula.
%        -- tol: Tolerancia  requerida.
%        -- N: Máximo número de iteraciones.

% OUTPUT: -- w: aproximaciones para u(x,y) que es la solución
%              de la ecuación diferencial de segundo orden.
%        -- l: número de iteraciones requeridas para obtener 
%              la solución.

function [] = PEFD(f,g,a,b,c,d,m,n,tol,N)

% h y k representan el tamaño de paso para la construcción de los
% puntos x_i e y_j respectivamente.
h = (b-a)/n;
k = (d-c)/m;

% Se crean x & y que son vectores que contienen las coordenadas  
% para los distintos puntos de la cuadrìcula que utiliza 
% el método.
x = zeros(1,n-1);
for i = 1:n-1
    x(i) = a+i*h;
end
y = zeros(1,m-1);
for j = 1:m-1
    y(j) = c+j*k;
end

% w es la matriz de ceros que luego se irá actualizando con 
% las aproximaciones w_ji de la solución u(x_i,y_j).
w = zeros(n-1,m-1);

% Se calculan lambda y mu para hacer más eficiente el código en térmios de
% las operaciones que se realizan. 
lambda = (h^2)/(k^2);
mu = 2*(1+lambda);

% El bucle while se utiliza para realizar iteraciones del método
% de Gauss-Seidel Optimizado
l = 1;
while l <= N     
    % Aquí se actualiza el elemento 1 de la última columna de w
    z = (-(h^2)*f(x(1),y(m-1))+g(a,y(m-1))+lambda*g(x(1),d)+lambda*w(1,m-2)+w(2,m-1))/mu;
    NORM = norm(z-w(1,m-1),inf);
    w(1,m-1)=z;
    
    % El siguiente for actualiza los valores de w, desde el segundo
    % al penúltimo elemento de w, de la última columna.
    for i = 2:n-2
        z = (-(h^2)*f(x(i),y(m-1))+lambda*g(x(i),d)+w(i-1,m-1)+w(i+1,m-1)+lambda*w(i,m-2))/mu;
        if norm(w(i,m-1)-z,inf) > NORM
            NORM = norm(w(i,m-1)-z,inf);
        end
        w(i,m-1)=z;
    end
    
    % Aquí se actualiza el último valor de la última columna de w.
    z = (-(h^2)*f(x(n-1),y(m-1))+g(b,y(m-1))+lambda*g(x(n-1),d)+w(n-2,m-1)+lambda*w(n-1,m-2))/mu;
    if norm(w(n-1,m-1)-z,inf) > NORM
        NORM = norm(w(n-1,m-1)-z,inf);
    end
    w(n-1,m-1) = z;
    % El siguiente for actualiza los valores de w, desde el segundo 
    % al penúltimo, de la todas las fila.
    for j = 2:m-2
        % Aquí actualizan los valores de w, desde el segundo al penúltimo, 
        % de la primera fila
        z = (-(h^2)*f(x(1),y(j))+g(a,y(j))+lambda*w(1,j+1)+lambda*w(1,j-1)+w(2,j))/mu;
        if norm(w(1,j)-z,inf) > NORM
            NORM = norm(w(1,j)-z,inf);
        end
        w(1,j) = z;
        
        % Este for actualiza los valores de w, desde el segundo al penúltimo, 
        % de las filas 2 a la n-2
        for i = 2:n-2
            z = (-(h^2)*f(x(i),y(j))+w(i-1,j)+lambda*w(i,j+1)+w(i+1,j)+lambda*w(i,j-1))/mu;
            if norm(w(i,j)-z,inf) > NORM
                NORM = norm(w(i,j)-z,inf);
            end
            w(i,j) =z;
        end
        % Aquí se actualizan los valores de w, desde el segundo al penúltimo,
        % de la última fila
        z = (-(h^2)*f(x(n-1),y(j))+g(b,y(j))+w(n-2,j)+lambda*w(n-1,j+1)+lambda*w(n-1,j-1))/mu;
        if norm(w(n-1,j)-z,inf) > NORM
            NORM = norm(w(n-1,j)-z,inf);
        end
        w(n-1,j) =z;
    end
    
    % Aquí se actualiza el primer valor de la primera fila y primera columna de w. fila 1 y columna 1 de w.
    z = (-(h^2)*f(x(1),y(1))+g(a,y(1))+lambda*g(x(1),c)+lambda*w(1,2)+w(2,1))/mu;
    if norm(w(1,1)-z,inf) > NORM
        NORM = norm(w(1,1)-z,inf);
    end
    w(1,1) =z;
    
    % El siguiente for actualiza los valores de w, desde el segundo al penúltimo
    % de la primera columna.
    for i = 2:n-2
        z = (-(h^2)*f(x(i),y(1))+lambda*g(x(i),c)+w(i-1,1)+lambda*w(i,2)+w(i+1,1))/mu;
        if norm(w(i,1)-z,inf) > NORM
            NORM = norm(w(i,1)-z,inf);
        end    
        w(i,1) =z;
    end
    
    % Aquí se actualiza el último valor de la ultima columna de la matriz w.
    z = (-(h^2)*f(x(n-1),y(1))+g(b,y(1))+lambda*g(x(n-1),c)+w(n-2,1)+lambda*w(n-1,2))/mu;
    if norm(w(n-1,1)-z,inf) > NORM
        NORM = norm(w(n-1,1)-z,inf);
    end
    w(n-1,1) =z;
    
    % Este if verifica que se ha alcanzado la tolerancia deseada.
    if NORM <= tol
        break
    end
    
    % Aquí se actualiza el contador del bucle while.
    l = l+1;
end

% Este último if envía un mensaje de error si se ha revasado el número
% de iteraciones deseado.
if l > N
    error('El numero maximo de iteraciones fue excedido') 
end
I=zeros(m-1,1);
J=zeros(m-1,1);
%imprime los títulos de cada columna de la tabla resultante
fprintf('      i         j         xi        xj      w(i,j) \n')
for j=1:m-1 %se llena los valores para los índices j en la tabla resultante
    J(j)=j;
end
for i=1:n-1 %se llenan los vectores con las coordenadas de x e y
    x1=zeros(m-1,1);
    for j=1:m-1
        I(j)=i; %se llena los valores para los indices i en la tabla resultante
        x1(j)=x(i);
    end
    %imprime los elemento de la tabla
    disp([ I  J  x1 transpose(y) transpose(w(i,:)) ])
end
% muestra el número de iteraciones
fprintf('\n Y el numero de iteraciones necesarias fue: %d ', l) 
end % termina el codigo

% Ejecución del ejemplo 2 , de la sección 12.1 del libro de Burden
% f = @(x,y)x*exp(y);
% g = @(x,y)x*exp(y);
% a=0; b=2; c=0; d=1; m=5; n=6; tol=1e-10; N=100; 
% PEFD(f,g,a,b,c,d,m,n,tol,N)

