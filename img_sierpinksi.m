

N = 10000;
ds_factor = 14;
Nds = round(N/ds_factor);

x = ones(1,N);
y = ones(1,N);
t = zeros(Nds);

for i = 2:N
    k = ceil(rand*3);
    x(i) = x(i-1)/2 + (k-1)*.25;
    y(i) = y(i-1)/2 + (k==2)*.5;
end

figure(3), clf
plot(x,y,'.')


figure(4), clf
t(sub2ind([Nds,Nds],ceil(Nds*y),ceil(Nds*x))) = 1;
imagesc(t), colormap gray, axis xy


