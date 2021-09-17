function [h]=altvertido(H1, Q, g, b, z, p)
  
  h=ones(1,numel(Q));;
  hant=zeros(1,numel(Q));;
  
  while max(abs(h-hant)) > 1e-6 
    hant = h;
    h= H1-(Q./(b.*(hant+p)+z.*(h+p).^2)).^2./(2*g);
  endwhile
  
endfunction