% Aforador de pared gruesa en canal de sección trapecial

close all;
clear, clc;

%{
y: calado
b: ancho de solera
z: talud de los quijeros
w=b·y+z·y²: sección de la corriente
p=b+2y·sqrt(1+z²): perímetro mojado de la corriente
B=b+2z·y: ancho superficial de la corriente
%}

addpath('src') %Ruta a la carpeta de archivos con las funciones

% Datos (expresados en el SI)
Q = 0:0.1:5; %Vector de valores de caudal
g = 9.80665; %Gravedad
bA = 3.44; %Ancho de solera en la aproximación al aforador
zAi = 0.2427; %Talud quijero (o cajero) izquierdo en la aproximación
zAd = 0.2064; %Talud quijero (o cajero) derecho en la aproximación
p = 1.1165; %Altura del realce de solera
I0B = 0; %Pendiente sobre el tramo horizontal del aforador
bB = bA+p*(zAi+zAd);%3.96; %Ancho solera sobre aforador --> ¡Ojo!, puede haber alguna discrepancia
zBi = zAi; %Talud quijero (o cajero) izquierdo sobre el aforador
zBd = zAd; %Talud quijero (o cajero) derecho sobre el aforador
L = 3.2; %Longitud del aforador
n = 0.012; %Coeficiente n de Manning
I0 = 0.016; %Pendiente topográfica de la rasante del canal aguas abajo del aforador
Lmodular = 0.75; %Valor de la relación entre la energía aguas abajo y la mínima sobre el resalto, en relación a la rasante del aforador

[yc, H0min] = ycritico(Q, g, bB, zBi, zBd);

%Cálculo de las pérdidas de carga. Seleccionar una opción, a) o b)

%{
% a) Rozamiento supuesto tramo uniforme con velocidad crítica
hf=L.*(Q.*n.*(bB+yc.*(sqrt(1+zBi)+sqrt(1+zBd))).^(2/3)./(bB.*yc+(zBi+zBd)/2.*yc.^2).^(5/3)).^2;
%}

% b) Rozamiento con Régimen Gradualmente Variado
[y1p, H01p]=RGV(Q, g, I0B, bB, zBi, zBd, L, n, yc, H0min);%Tramo horizontal sobre el aforador
hf=H01p-H0min;

H1 = H0min + hf; %Ecuación de la energía. Podrían añadirse pérdidas de carga H1 = H0min + hf

[h] = altvertido(H1, Q, g, bA, zAi, zAd, p);

Resultado=[Q; h; yc; H0min; h+p];

% Ajuste mediante mínimos cuadrados
% [p,s,mu]=polyfit(x,y,n)
[p10,s10]=polyfit(log(Resultado(2,2:end)), log(Resultado(1,2:end)), [true true]);

errorp10=transpose(sqrt(diag(s10.C)/s10.df)*s10.normr);
err_rel_p10=errorp10./p10;

% Se muestran los coeficientes de ajuste y sus errores de determinación en la consola
disp('  c1  log(c0) ')
disp(p10)
disp('  +-  +-  +-  ')
disp(errorp10)
disp('Error relativo')
disp(err_rel_p10)
disp(' c0 ')
disp(exp(p10(2)))
disp('  +-  ')
disp(exp(p10(2)+errorp10(2))-exp(p10(2)))
disp('--------')

Q_est=exp(p10(1,2)).*Resultado(2,:).^p10(1,1);

figure()
plot(Resultado(1,:),Resultado(5,:),'*')
titulo = sprintf('Ajuste potencial Q = %d·h^{%d}', exp(p10(2)),p10(1));
title(titulo)
hold on
plot(Q_est, Resultado(5,:))
xlabel('Q(m^{3}/s)')
%ylabel('h(m) = y-p(m)')
ylabel('y(m) = h+p(m)')
hold off

yNormal=caladoNormal(Q, g, I0, bA, zAi, zAd, n);

figure()
title('Límite modular según n')
xlabel('Q(m^{3}/s)')
ylabel('L')
ylim([0 1.2])
hold on
for n=0.012:0.002:0.03
  yNormal=caladoNormal(Q, g, I0, bA, zAi, zAd, n);
  H0=yNormal+(Q./(bA*yNormal+(zAi+zAd)/2*yNormal.^2)).^2/2/g;
  L=(H0-p)./Resultado(4,:);
  plot(Q, L)
  etiqueta=sprintf('n=%d', n);
  text(Q(end), L(end), etiqueta)
endfor
plot([Q(1) Q(end)], [Lmodular Lmodular])
hold off
