function [y1p, H01p]=RGV(Q, g, I0, b, zi, zd, L, n, yc, H0min)
  
  Ic=(n*Q.*(b+(sqrt(1+zi)+sqrt(1+zd))*yc).^(2/3)./(b*yc+(zi+zd)/2*yc.^2).^(5/3)).^2;
  
  Dy1p=0.001*ones(1,numel(Q)); % Delta de y inicial 1 mm
  Dx=zeros(1,numel(Q));
    
  while max(abs(L-Dx)) > 1e-6
            
    y1p=yc+Dy1p;
    H01p=y1p+(Q./(b.*y1p+(zi+zd)/2.*y1p.^2)).^2./2./g;
    I1p=(n*Q.*(b+y1p.*(sqrt(1+zi)+sqrt(1+zd))).^(2/3)./(b.*y1p+(zi+zd)/2.*y1p.^2).^(5/3)).^2;
    Dx=(H01p-H0min)./((Ic+I1p)./2-I0);
    
    correc = L./Dx;
    Dy1p = Dy1p.*correc.^(1/2);
    disp(mean(correc(2:end)))
    
  endwhile
    
endfunction