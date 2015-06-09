% plot lists of line definitions, like the ones
% returned from subCannyEdgeMEX()
function plotEdges(e,sym)
  if (nargin < 2), sym = 'b'; end
  n = length(e);
  for k=1:n
    pp = e(k).p;  % for now, it only returns points lists.
	  % someday we may add the direction of the line at each point
    j = mod(k-1,length(sym))+1;
    plot(pp(2,:),-pp(1,:),sym(j));
    if (k == 1), hold on; end
  end
end
