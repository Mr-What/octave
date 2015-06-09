% print out edges as paths in an SVG image
% returned from subCannyEdgeMEX()
function edge2svg(e,fNam)
  fid = fopen(fNam,'wt');
  fprintf(fid,'<?xml version="1.0" encoding="utf-8"  standalone="no"?>\n');
  fprintf(fid,'<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n');
  ax = getAxisExtent(e)
  dx = ax(2)-ax(1);
  dy = ax(4)-ax(3);
  x0 = floor(ax(1));  x1 = ceil(ax(2));
  y0 = floor(ax(3));  y1 = ceil(ax(4));
  dx = x1-x0+1;
  dy = y1-y0+1;
  fprintf(fid,'<svg width="%d" height="%d"\n\tviewBox="%d %d %d %d"\n',...
	  dx,dy,x0,y0,dx,dy);
  fprintf(fid,'\tpreserveAspectRatio="xMinYMin meet"\n');
  fprintf(fid,'\txmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">\n');

  n = length(e);
  for k=1:n
    pp = e(k).p;
    m = size(pp,2);
    fprintf(fid,'<path style="stroke:#0000dd; fill:none" d="M%.2f,%.2f\n',pp(2,1),pp(1,1));
    for j=2:m-1
      fprintf(fid,'L%.2f,%.2f\n',pp(2,j),pp(1,j));
    end
    d = pp(:,end) - pp(:,1);
    d = norm(d)
    if (d > 2)
      fprintf(fid,'L%.2f,%.2f" />\n',pp([2 1],end));
    else
      fprintf(fid,'Z" />\n');
    end
  end

  fprintf(fid,'</svg>\n');
  fclose(fid);
end

function ax = getAxisExtent(e)
  n = length(e);
  xRange = range(e(1).p(2,:));
  yRange = range(e(1).p(1,:));
  for k=2:n
    xRange = range([xRange,e(k).p(2,:)]);
    yRange = range([yRange,e(k).p(1,:)]);
  end
  ax = [xRange,yRange]
end

