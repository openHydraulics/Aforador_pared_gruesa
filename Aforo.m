# Aforador de pared gruesa en canal de sección trapecial
#{
y: calado
b: ancho de solera
z: talud de los quijeros
w=b·y+z·y²: sección de la corriente
p=b+2y·sqrt(1+z²): perímetro mojado de la corriente
B=b+2z·y: ancho superficial de la corriente
#}

addpath('src') %Ruta a la carpeta de archivos con las funciones

% Datos (expresados en el SI)
Q = 0:0.1:5; %Caudal
g = 9.80665; %Gravedad
bB = 1; %Ancho solera sobre aforador
zB = 0.58; %Talud quijeros (o cajeros) sobre el aforador
p = 0.2; %Altura del realce de solera
bA = 2; %Ancho de solera en la aproximación al aforador
zA = 0.58; %Talud de quijeros en la aproximación

[yc, H0min] = ycritico(Q, g, bB, zB);

H1 = H0min; %Ecuación de la energía. Podrían añadirse pérdidas de carga H1 = H0min + hf

[h] = altvertido(H1, Q, g, bA, zA, p);

Resultado=[Q; h; yc; H0min; h+p];

# Ajuste mediante mínimos cuadrados
% [p,s,mu]=polyfit(x,y,n)
[p10,s10]=polyfit(log(Resultado(2,2:end)), log(Resultado(1,2:end)), [true true]);

errorp10=transpose(sqrt(diag(s10.C)/s10.df)*s10.normr);
err_rel_p10=errorp10./p10;

# Se muestran los coeficientes de ajuste y sus errores de determinación en la consola
disp('  c1  log(c0) ')
disp(p10)
disp('  +-  +-  +-  ')
disp(errorp10)
disp('Error relativo')
disp(err_rel_p10)
disp('--------')

Q_est=exp(p10(1,2)).*Resultado(2,:).^p10(1,1);

# Representación de los cálculos
figura=figure();
plot(Resultado(1,:),Resultado(2,:),'*')
xlabel('Q(m^3/s)')
ylabel('h(m) = y-p(m)')
title('Ajuste potencial Q = c0·h^c^1')
hold on
plot(Q_est,Resultado(2,:))

# Guarda la figura como jpg con el nombre figura.jpg
%print(gcf,'figura.jpg', '-djpg','-r0')
%print(figura,'figura.jpg', '-djpg','-r0')