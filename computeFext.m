function [Fext] = computeFext(n_dof,n_el,n_ne,n_i,Td,Me,g,Fel,x,L1)
% Compute engine's weight (point force)

Fext=zeros(n_dof,1);
for i=1:n_el
    if (x(i)==L1)
        Fext(2*i-1)=-Me*g;
    else 
        if(x(i)>L1) && (x(i-1)<L1)
        Fext(2*(i-1)-1)=-Me*g/2;
        Fext(2*i-1)=-Me*g/2;
        end
    end
end

% This was in assemblyKG before (unnecessary)
for e=1:n_el
    for i=1:(n_ne*n_i)
        I=Td(e,i);
        Fext(I)=Fext(I)+Fel(i,e);
    end
end

end