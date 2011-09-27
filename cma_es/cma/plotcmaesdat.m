function plotcmaesdat(figNb, filenameprefix, filenameextension, objectvarname)
% PLOTCMAESDAT;
% PLOTCMAES(FIGURENUMBER_iBEGIN_iEND, FILENAMEPREFIX, FILENAMEEXTENSION, OBJECTVARNAME);
%   plots output from CMA-ES, e.g. cmaes.m, Java class CMAEvolutionStrategy...
%   mod(figNb,100)==1 plots versus iterations.
%
% PLOTCMAES([101 300]) plots versus iteration, from iteration 300.
% PLOTCMAES([100 150 800]) plots versus function evaluations, between iteration 150 and 800.
%
% Upper left subplot: blue/red: function value of the best solution in the
%   recent population, cyan: same function value minus best
%   ever seen function value, green: sigma, red: ratio between
%   longest and shortest principle axis length which is equivalent
%   to sqrt(cond(C)).
% Upper right plot: time evolution of the distribution mean (default) or
%   the recent best solution vector.
% Lower left: principle axes lengths of the distribution ellipsoid,
%   equivalent with the sqrt(eig(C)) square root eigenvalues of C.
% Lower right: magenta: minimal and maximal "true" standard deviation
%   (with sigma included) in the coordinates, other colors: sqrt(diag(C))
%   of all diagonal elements of C, if C is diagonal they equal to the
%   lower left.
%
% Files [FILENAMEPREFIX name FILENAMEEXTENSION] are used, where
%   name = axlen, OBJECTVARNAME (xmean|xrecentbest), fit, or stddev.
%
manual_mode = 1;

  if nargin < 1 || isempty(figNb)
    figNb = 325;
  end
  if nargin < 2 || isempty(filenameprefix)
    filenameprefix = 'outcmaes';
  end
  if nargin < 3 || isempty(filenameextension)
    filenameextension = '.dat';
  end
  if nargin < 4 || isempty(objectvarname)
    objectvarname = 'xmean';
    objectvarname = 'xrecentbest';
  end

  % load data
  d.x = load([filenameprefix objectvarname filenameextension]);
  % d.x = load([filenameprefix 'xmean' filenameextension]);
  % d.x = load([filenameprefix 'xrecentbest' filenameextension]);
  d.f = load([filenameprefix 'fit' filenameextension]);
  d.std = load([filenameprefix 'stddev' filenameextension]);
  d.D = load([filenameprefix 'axlen' filenameextension]);

  % interpret entries in figNb for cutting out some data
  if length(figNb) > 1
    iend = inf;
    istart = figNb(2);
    if length(figNb) > 2
      iend = figNb(3);
    end
    figNb = figNb(1);
    d.x = d.x(d.x(:,1) >= istart & d.x(:,1) <= iend, :);
    d.f = d.f(d.f(:,1) >= istart & d.f(:,1) <= iend, :);
    d.std = d.std(d.std(:,1) >= istart & d.std(:,1) <= iend, :);
    d.D = d.D(d.D(:,1) >= istart & d.D(:,1) <= iend, :);
  end

  % set up figure window
  if manual_mode
    figure(figNb);  % just create and raise the figure window
  else  % try to prevent raise
    if evalin('caller', 'countiter') <= 2 && evalin('caller', 'irun') == 1
      figure(324);
    elseif gcf ~= 324
      if ismember(324, findobj('Type', 'figure'))
        set(0, 'CurrentFigure', 324);  % prevents raise of figure window
      else
        figure(324);
      end
    end
  end

  % decide for x-axis
  iabscissa = 2; % 1== versus iterations, 2==versus fevals
  if mod(figNb,100) == 1
    iabscissa = 1; % a short hack
  end
  if iabscissa == 1
    xlab ='iterations';
  elseif iabscissa == 2
    xlab = 'function evaluations';
  end

  if size(d.x, 2) < 1000
    minxend = 1.03*d.x(end, iabscissa);
  else
    minxend = 0;
  end

  % plot fitness etc
  foffset = 1e-99;
  dfit = d.f(:,6)-min(d.f(:,6));
  [ignore idxbest] = min(dfit);
  dfit(dfit<1e-98) = NaN;
  subplot(2,2,1); hold off;
  dd = abs(d.f(:,7:8)) + foffset;
  dd(d.f(:,7:8)==0) = NaN;
  semilogy(d.f(:,iabscissa), dd, '-k'); hold on;
  % additional fitness data, for example constraints values
  if size(d.f,2) > 8
    dd = abs(d.f(:,9:end)) + 10*foffset;  % a hack
    % dd(d.f(:,9:end)==0) = NaN;
    semilogy(d.f(:,iabscissa), dd, '-m'); hold on;
    if size(d.f,2) > 12
      semilogy(d.f(:,iabscissa),abs(d.f(:,[7 8 11 13]))+foffset,'-k'); hold on;
    end
  end

  idx = find(d.f(:,6)>1e-98);  % positive values
  if ~isempty(idx)  % otherwise non-log plot gets hold
    semilogy(d.f(idx,iabscissa), d.f(idx,6)+foffset, '.b'); hold on;
  end
  idx = find(d.f(:,6) < -1e-98);  % negative values
  if ~isempty(idx)
    semilogy(d.f(idx, iabscissa), abs(d.f(idx,6))+foffset,'.r'); hold on;
  end
  semilogy(d.f(:,iabscissa),abs(d.f(:,6))+foffset,'-b'); hold on;
  semilogy(d.f(:,iabscissa),dfit,'-c'); hold on;
  semilogy(d.f(:,iabscissa),(d.f(:,4)),'-r'); hold on; % AR
  semilogy(d.std(:,iabscissa), [max(d.std(:,6:end), [], 2) ...
                      min(d.std(:,6:end), [], 2)], '-m', 'linewidth', 2); % max,min std
  [maxval, imax] = max(d.std(end,6:end));
  [minval, imin] = min(d.std(end,6:end));
  text(d.std(end,iabscissa), maxval, sprintf('%d:%.0e', imax, maxval));
  text(d.std(end,iabscissa), minval, sprintf('%d:%.0e', imin, minval));

  semilogy(d.std(:,iabscissa),(d.std(:,3)),'-g'); % sigma
  % plot best f
  semilogy(d.f(idxbest,iabscissa),min(dfit),'*c'); hold on;
  semilogy(d.f(idxbest,iabscissa),abs(d.f(idxbest,6))+foffset,'*r'); hold on;

  ax = axis;
  ax(2) = max(minxend, ax(2));
  axis(ax);

  yannote = 10^(log10(ax(3)) + 0.05*(log10(ax(4))-log10(ax(3))));
  text(ax(1), yannote, ...
       [ 'f=' num2str(d.f(end,6), '%.15g') ]);

  title('blue:abs(f), cyan:f-min(f), green:sigma, red:axis ratio');
  grid on;
%  ax([3,4]) = [1e-9, 1e5];
%  axis(ax);

  % plot x-values
  subplot(2,2,2); hold off;
  plot(d.x(:,iabscissa), d.x(:,6:end),'-'); hold on;
  ax = axis;
  ax(2) = max(minxend, ax(2));
  axis(ax);

  % add some annotation lines
  [ignore idx] = sort(d.x(end,6:end));
  % choose no more than 25 indices
  idxs = round(linspace(1, size(d.x,2)-5, min(size(d.x,2)-5, 25)));
  yy = repmat(NaN, 2, size(d.x,2)-5);
  yy(1,:) = d.x(end, 6:end);
  yy(2,idx(idxs)) = linspace(ax(3), ax(4), length(idxs));
  plot([d.x(end,iabscissa) ax(2)], yy, '-');
  plot(repmat(d.x(end,iabscissa),2), [ax(3) ax(4)], 'k-');
  for i = idx(idxs)
    text(ax(2), yy(2,i), ...
         ['x(' num2str(i) ')=' num2str(yy(1,i), '%.3g')]);
  end
  title(['Object Variables, ' objectvarname(2:min(end,7)) ' (' num2str(size(d.x, 2)-5) '-D)']);grid on;

  subplot(2,2,3); hold off; semilogy(d.D(:,iabscissa), d.D(:,6:end), '-');
  ax = axis;
  ax(2) = max(minxend, ax(2));
  axis(ax);
  title('Principle Axes Lengths');grid on;
  xlabel(xlab);

  subplot(2,2,4); hold off;
  % semilogy(d.std(:,iabscissa), d.std(:,6:end), 'k-'); hold on;
  % remove sigma from stds
  d.std(:,6:end) = d.std(:,6:end) ./ (d.std(:,3) * ones(1,size(d.std,2)-5));
  semilogy(d.std(:,iabscissa), d.std(:,6:end), '-'); hold on;
  if 11 < 3  % max and min std
    semilogy(d.std(:,iabscissa), [d.std(:,3).*max(d.std(:,6:end), [], 2) ...
                        d.std(:,3).*min(d.std(:,6:end), [], 2)], '-m', 'linewidth', 2);
    maxval = max(d.std(end,6:end));
    minval = min(d.std(end,6:end));
    text(d.std(end,iabscissa), d.std(end,3)*maxval, sprintf('max=%.0e', maxval));
    text(d.std(end,iabscissa), d.std(end,3)*minval, sprintf('min=%.0e', minval));
  end
  ax = axis;
  ax(2) = max(minxend, ax(2));
  axis(ax);
  % add some annotation lines
  [ignore idx] = sort(d.std(end,6:end));
  % choose no more than 25 indices
  idxs = round(linspace(1, size(d.x,2)-5, min(size(d.x,2)-5, 25)));
  yy = repmat(NaN, 2, size(d.std,2)-5);
  yy(1,:) = d.std(end, 6:end);
  yy(2,idx(idxs)) = logspace(log10(ax(3)), log10(ax(4)), length(idxs));
  semilogy([d.std(end,iabscissa) ax(2)], yy, '-');
  semilogy(repmat(d.std(end,iabscissa),2), [ax(3) ax(4)], 'k-');
  for i = idx(idxs)
    text(ax(2), yy(2,i), [' ' num2str(i)]);
  end
  title('Standard Deviations in Coordinates divided by sigma');grid on;
  xlabel(xlab);

  if figNb ~= 324
    % zoom on;  % does not work in Octave
  end
  drawnow;
