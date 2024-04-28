% Test 1D geometries
% Task: Create grid to discretize [0,1] by using different calls,
% and plot a function.

% A function for testing plot
f = @(x) sin(2*pi*x)+sin(4*pi*x);

% All valid calls
I = Interval([0 1]);
figure
I.plot(f(I.x),'Marker','.');

I = Interval([0 1],0.1);
figure
I.plot(f(I.x),'Marker','.');


I = Interval([0 0.2 0.5 1]);
I.refineMesh;figure
I.plot(f(I.x),'Marker','.');

I = Interval([0 0.2 0.5 1],0.0625);
figure
I.plot(f(I.x),'Marker','.');

pause(5)
close all

% Force some typical mistakes that user maybe will do.
% Empty constructor not allowed
try
  I = Interval();
catch ME
  fprintf(['Error: ',ME.message,'\n'])
end
% Missformed interval
try
  I = Interval([1 1]);
catch ME
  fprintf(['Error: ',ME.message,'\n'])
end

% Matrix instead vector
try
  I = Interval([0 1; 2 3]);
catch ME
  fprintf(['Error: ',ME.message,'\n'])
end
% Scalar [18] insted of [1 8]
try
  I = Interval([18]);
catch ME
  fprintf(['Error: ',ME.message,'\n'])
end
% Char 'H' instead of H.
try
  I = Interval([0 1],'H');
catch ME
  fprintf(['Error: ',ME.message,'\n'])
end
% Negative  meshwidth
try
  I = Interval([0 1],-1);
catch ME
  fprintf(['Error: ',ME.message,'\n'])
end
% Too many arguments 1
try
  I = Interval([0 1],0.1,1);
catch ME
  fprintf(['Error: ',ME.message,'\n'])
end
% Too many arguments 2: Wrote 0,1 instead  0:1
try
  I = Interval(0,1,0.1);
catch ME
  fprintf(['Error: ',ME.message,'\n'])
end
