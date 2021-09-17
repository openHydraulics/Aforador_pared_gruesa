function [yc, H0min]=ycritico(Q, g, b, z)
  
  yc=ones(1,numel(Q));
  ycant=zeros(1,numel(Q));
  
  while max(abs(yc-ycant)) > 1e-6 
    ycant = yc;
    yc= ((Q.^2.*(b+2*z.*ycant))./(g*(b+z*ycant).^3)).^(1/3);
  endwhile
  
  H0min = yc+(b.*yc+z.*yc.^2)/2./(b+2*z.*yc);
  
endfunction
