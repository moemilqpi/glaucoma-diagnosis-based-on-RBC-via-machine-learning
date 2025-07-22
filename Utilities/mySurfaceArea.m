function surface = mySurfaceArea(X,Y,H,pixel)
    dx = pixel;
    dy = pixel;
    dA = zeros(length(X),1);
    [Fx,Fy] = gradient(H,dx,dy);
    for i=1:length(X)
        dA(i) = dx*dy*sqrt(1+(Fx(X(i),Y(i)))^2+(Fy(X(i),Y(i)))^2);
    end
    surface = sum(dA);
end