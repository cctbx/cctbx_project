//_______________________________________________________
//_______________________________________________________
//
function res = minnan(in)

  res = min(in(~isnan(in)));

endfunction

//_______________________________________________________
//_______________________________________________________
//
function plotcmaesdat(name_prefix, fignb, name_extension, object_variables_name)
//
// plots data from CMA-ES, so far for Scilab and Java.
// defaults: name_prefix='outcmaes', name_extension='.xml'
//           fignb=326

// Can be called by Plot Now button, therefor the strange
// arrangement of input elements.

defaults.name_prefix = 'outcmaes';
defaults.name_extension = '.dat';
defaults.object_variables = 'xmean'; // or 'xrecentbest'

if ~isdef('name_prefix', 'local')
  name_prefix = defaults.name_prefix;
  if isdef('filenameprefix') // inherit from upper environment
    name_prefix = filenameprefix;
  end
elseif isempty(name_prefix)
  name_prefix = defaults.name_prefix;
elseif name_prefix == 1 // called from "Plot Now" button
  name_prefix = get(gcf(), "figure_name");
end
if ~isdef('name_extension', 'local')
  name_extension = defaults.name_extension;
elseif isempty(name_extension)
  name_extension = defaults.name_extension;
end
if ~isdef('fignb', 'local')
  fignb = 326;
end

if ~isdef('object_variables_name', 'local')
  object_variables_name = defaults.object_variables;
end


  names = ['axlen', 'fit', 'stddev', object_variables_name];
//  names = ['axlen', 'fit', 'stddev', 'xmean', 'xrecentbest'];
  data = tlist(['data', names]);
// columns are generation feval something y-values

  // read data from files into tlist data
  for name = names
    [fid, err] = mopen([name_prefix + name + name_extension]);
    if err ~= 0
      warning('File ' + [name_prefix + name + name_extension] + ...
          ' could not be opened');
      data(name) = [];
    else
      if 1 < 3  // reading quick and dirty
        mclose(fid);
        data(name) = fscanfMat([name_prefix + name + name_extension]);
      else // this can become a more fail save version of reading the data
        // TODO: check whether this also works if several DATA entries are present
        s = mgetl(fid); // read complete file
        mclose(fid);
        // remove headings up to < DATA ... >
        idx = grep(s, 'DATA'); // TODO: make fail save by searching for < DATA * >
        idx2 = grep(s, '>'); // lines where > character is found
        idx2 = idx2(idx2>=idx(1));
        if ~isempty(idx2)
          s = s(idx2(1)+1:$); // TODO keep eventually also the header stuff
        end

        // remove trailing non numerical stuff
        idx = grep(s, '<');
        if ~isempty(idx)
          s = s(1:idx(1)-1);
        end
        data(name) = evstr(s);
        //       data(name) = mfscanf(-1, fid);
      end
    end
  end

// TODO the legend interferes with large negativ fitnesses

// comments:
// data.*(:,1) == iteration number
// data.*(:,2) == function evaluations
// data.*(:,6:$) == the vector data like xmean
// data.fit(:,6) == recent best fitness

  flg_raw_stds = %F; // remove sigma from stds
  if flg_raw_stds & ~isempty(data.stddev)
    for i = 6:size(data.stddev,2)
      data.stddev(:,i) = data.stddev(:,i) ./ data.stddev(:,3);
    end
  end

  scf(fignb); // set current figure
  drawlater;
  clf; // is needed only because of the annotation of variables


  //////////////////
  subplot(2,2,1);
  if ~isempty(data.fit)
// TODO remove this (it is for testing purpose only)
// data.fit(:,6) = data.fit(:,6) - 2e5;
    xgrid;
    fmin = min(data.fit(:,6));
    idxmin = find(data.fit(:,6) == fmin);

    xtitle("Function Value (fval, fval minus f_min), Sigma (g), Axis Ratio (r)", ...
        'f_recent='+sprintf('%25.18e', data.fit($, 6)), ...
        'log10(abs(value))');

    // plot legend for best function value
    plot(data.fit(:,2)(idxmin), log10(abs(fmin) + 1e-49), 'r*');
    legend('f_best=' + sprintf('%25.18e', fmin));
    data.fit(idxmin,6) = %nan;

    // plot abs function value in blue, red, and black
    idxpos = data.fit(:,6) > 0;
    idxneg = data.fit(:,6) < 0; // removes nan and zeros

    if or(idxpos) // if any entry in idxpos is true
      plot(data.fit(idxpos,2), ...
          log10(0e-49 + abs(data.fit(idxpos,6))), 'b.');
    end

    if or(idxneg)
      plot(data.fit(idxneg,2), ...
          log10(0e-49 + abs(data.fit(idxneg,6))), 'k.');
    end

    // median and worst fitness
    if size(data.fit,2) > 6
      plot(data.fit(:,2), log10(abs(data.fit(:,7:8)) + 1e-99), 'k-');
    end
    // min and max stddev
    if size(data.fit,2) > 12
      plot(data.fit(:,2), log10(abs(data.fit(:,[11 13]))), 'k-');
    end

    // plot function value differences disregarding all fmins
    // careful: unfortunately log10 cannot handle %nan (in some versions)
    idx = ~isnan(data.fit(:,6));
    plot(data.fit(idx,2), log10(data.fit(idx,6) - fmin + 1e-49), 'c-');

    // plot marker(s) for best function value _over_ the graph now, see legend above
    plot(data.fit(:,2)(idxmin), log10(abs(fmin) + 1e-49), 'r*');
    plot(data.fit(:,2)(idxmin), log10(minnan(data.fit(:,6)-fmin) + 1e-49)*ones(idxmin), 'r*');

  end // data.fit is not empty

  // plot sigma
  if ~isempty(data.stddev)
    plot(data.stddev(:,2), log10(data.stddev(:,3)), 'g');
  else
    plot(data.fit(:,2), log10(data.fit(:,3)), 'g');
  end

  // plot axis ratio
  plot(data.fit(:,2), log10(data.fit(:,4)), 'r');

  //////////////////
  subplot(2,2,2);
  name = object_variables_name; // 'xmean' or 'xrecentbest'
  if ~isempty(data(name))
    xgrid;
    xtitle('Object Variables (' + name + ', ' + ...
        string(size(data(name),2)-5) + '-D)');

    plot(data(name)(:,2), data(name)(:,6:$));

    // annotations of variables with numbers
    Ngeno = size(data(name), 2) - 5;
    if Ngeno < 100

      yrange = get(gca(), "y_ticks");
      yrange = yrange.locations([1 $]);
      yrange(1) = yrange(1) + 1e-6*(yrange(2) - yrange(1)); // prevents widening of range
      yrange(2) = yrange(2) - 1e-6*(yrange(2) - yrange(1));
      xrange = get(gca(), "x_ticks");
      xrange = xrange.locations([1 $]);
      xlast = max([1.07 * data(name)($,2) xrange(2)]);
      [sorted_x, idx] = sort(data(name)($,6:$));
      idx2 = [];
      idx2(idx) = (Ngeno-1:-1:0)';
      plot(data(name)($,2)*ones(1,2), yrange, 'k', 'linewidth', 1);
      plot([data(name)($,2); xlast]', ...
          [data(name)($, 6:$) ; ...
              yrange(1) + (idx2')*(yrange(2)-yrange(1))/(Ngeno-1)]);

      set(gca(), "clip_state", "off"); // numbers visible everywhere
      str = ' ';
      for i = 1:Ngeno
        if i > 9
          str = '';
        end
        xstring(xlast, yrange(1) + ...
            (yrange(2)-yrange(1))*(idx2(i)/(Ngeno-1) - 0.3/max(10,Ngeno)), ...
            [str string((i))]);
      end
    end // annotation

  end // isempty data. xmean

  //////////////////
  subplot(2,2,3);
  if ~isempty(data.axlen)
    xgrid;
    xtitle("Principle Axis Lengths", 'function evaluations', 'log10(value)');
    plot(data.axlen(:,2), log10(data.axlen(:,6:$)+1e-49));
  end

  //////////////////
  subplot(2,2,4);
  if ~isempty(data.stddev)
    xgrid;
    xtitle("Standard Deviations", 'function evaluations', 'log10(value)');

    plot(data.stddev(:,2), log10(data.stddev(:,6:$)));

    // annotations of variables with numbers
    Ngeno = size(data.stddev, 2) - 5;
    if Ngeno < 100
      yrange = get(gca(), "y_ticks");
      yrange = yrange.locations([1 $]);
      yrange(1) = yrange(1) + 1e-6*(yrange(2) - yrange(1));
      yrange(2) = yrange(2) - 1e-6*(yrange(2) - yrange(1));
      xrange = get(gca(), "x_ticks");
      xrange = xrange.locations([1 $]);
      xlast = max([1.07 * data.stddev($,2) xrange(2)]);
      [sorted_x, idx] = sort(data.stddev($,6:$));
      idx2 = [];
      idx2(idx) = (Ngeno-1:-1:0);
      plot(data.stddev($,2)*ones(1,2), yrange, 'k', 'linewidth', 1);
      plot([data.stddev($,2); xlast]', ...
          [log10(data.stddev($,6:$)); ...
              yrange(1) + (yrange(2)-yrange(1))*(idx2/(Ngeno-1))]);

      set(gca(), "clip_state", "off");
      str = ' ';
      idxvars = 1:Ngeno; // find(out.genopheno.scaling > 0);
      for i = 1:Ngeno
        if i > 9
          str = '';
        end
        xstring(xlast, yrange(1) + ...
            (yrange(2)-yrange(1)) * (idx2(i)/(Ngeno-1) - 0.3/max(Ngeno,10)), ...
            [str string(idxvars(i))]);
      end
    end
  end // ~isempty(data.stddev)

  drawnow;

endfunction
