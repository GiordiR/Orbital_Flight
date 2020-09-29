function df = optinc(x)
% x: theta1
f2=fopen('Rnorm.dat');
dat=fscanf(f2,'%g',inf); %legge la prima riga di Rnorm.dat e la scrive nel vettore Rc
fclose(f2)
R=dat(1:3,1);
di= dat(4,1);
df = (R(1)*sin(x))/sqrt(1+R(1)^2-2*R(1)*cos(x)) + ((R(2)^2)*R(3)*sin(di-x))/sqrt(R(2)^2+(R(1)^2)*(R(3)^2)-2*(R(2)^2)*R(3)*cos(di-x));
end

