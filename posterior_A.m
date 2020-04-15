function [val] = posterior_A(A,B,b,S,omega)


global n p T

A0 = construct_A(A);
%A0 = vec2mat(A',n);


  if det(A0)<Inf && isreal(det(A0)) && det(A0)>0
      val = -((T-p)*log(det(A0))-1/2*trace((S+((B-b)'/(omega))*(B-b))*(A0'*A0)));
  else
      val = 10000000;
  end


end


