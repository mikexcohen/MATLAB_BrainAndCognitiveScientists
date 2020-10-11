

% Book cover: symmetric matrix with images, different color on diagonal.



n=10;
i=-n:.1:n;
d=zeros(length(i));

for ii=1:length(i)
    for jj=1:length(i)
        d(ii,jj) = log( complex(i(ii),i(jj)) );
    end
end

d(~isfinite(d))=min(d(:));

surf(abs(d)+real(d)+imag(d)), axis xy, axis square, axis off
shading interp, rotate3d on, set(gcf,'color','k')

cmap=(1+[cos(linspace(0,pi*2,100)); sin(linspace(0,pi*2,100)); cos(linspace(0,pi*2,100))])/2;
% colormap(cmap')

circmap = (1+[cos(linspace(0,pi*2,100)); sin(linspace(0,pi*2,100)); cos(linspace(0,pi*2,100) + pi)])/2;
colormap(circmap')



view([ -125  20 ])
view([  -40 -50 ])
view([  -30  40 ])

% diagonal
view([ -15  50 ])
view([ -45  65 ])
view([ -160 60 ])


