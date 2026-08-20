function L = distancept(p1,p2)


D=p1-p2;
L = sqrt(D(:,1).^2+D(:,2).^2+D(:,3).^2);