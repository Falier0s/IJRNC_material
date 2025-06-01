function [T, y]=RungeKutta(Fun,t_span,y0,step_size,u)

T=t_span(1):step_size:t_span(2);
y = zeros(length(y0),length(T)); 
if size(u,2)==1
    u=u.*ones(1,size(T,2));
end

y(:,1)=y0;
    for i=1:(length(T)-1)                              % calculation loop
        k_1 = Fun(T(i),y(:,i),u(:,i));
        k_2 = Fun(T(i)+0.5*step_size,y(:,i)+0.5*step_size*k_1,u(:,i));
        k_3 = Fun((T(i)+0.5*step_size),(y(:,i)+0.5*step_size*k_2),u(:,i));
        k_4 = Fun((T(i)+step_size),(y(:,i)+k_3*step_size),u(:,i));
        y(:,i+1) = y(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*step_size;  % main equation
      
    end

end