function [sse,p] = fit2segLinear(bpoint,x,y)

% sort according to x-axis
[~,i]=sort(x);
x=x(i); y=y(i);

% make sure it's column-format
x=x(:); y=y(:);

bpoint = abs(round(bpoint));

% force breakpoint to stay in the range of the data
bpoint = min(bpoint,length(x));
bpoint = max(1,bpoint);

%% fit linear models

% first piece
x1 = [ones(bpoint,1) x(1:bpoint)];
y1 = y(1:bpoint);
b1 = (x1'*x1)\(x1'*y1);

% second piece
x2 = [ones(length(x)-bpoint,1) x(bpoint+1:end)];
y2 = y(bpoint+1:end);
b2 = (x2'*x2)\(x2'*y2);

%% evaluate model

% predicted data
yHat = [b1(1)+b1(2)*x1(:,2); b2(1)+b2(2)*x2(:,2)];

% predictors and sse (minimization objective)
sse = sum( (yHat-y).^2 ) / sum( y.^2 );
p = [ b1' b2' x(bpoint) sse ];

% optional plotting
plot(x,y,'ro',x,yHat,'k-'); drawnow; pause(.01)

%%
        