function y0=caladoNormal(Q, g, I0, b, zi, zd, n)
  
  y0=1.5*ones(1,numel(Q));
  y0ant=y0./2;
    
  while max(abs(y0-y0ant)) > 1e-6
    y0ant=y0;
    y0=(n*Q.*(b+y0ant.*(sqrt(1+zi)+sqrt(1+zd))).^(2/3)./(b+(zi+zd)/2.*y0ant).^(5/3)/I0).^(3/5);    
  endwhile
    
endfunction