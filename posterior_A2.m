function [val] = posterior_A2(A,B,b,S,omega)


global n p T

%A0 = vec2mat(A,n);
A0 = construct_A(A);
  if det(A0)<Inf && isreal(det(A0)) && det(A0)>0
      val = ((T-p)*log(det(A0))-1/2*trace((S+((B-b)'/(omega))*(B-b))*(A0'*A0)));
      val = (val);
  else
      val = 0;
  end

end


