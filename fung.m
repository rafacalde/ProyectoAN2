function [g] = fung(x,y)
%%%%%% ACA SE GENERAN LAS CONDICIONES DE FRONTERA
if x==0;
    g=0;
end
if x==2;
    g=2*exp(y);
end
if y==0;
    g=x;
end
if y==1;
    g=exp(1)*x;
end
end
