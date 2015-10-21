% Given an edge, like from subCannyEdgeMEX()
% check to see if the last few samples seem to overlap the first one(s)
% as in a loop.  If so, trim the last samples, and close the loop
% to the initial sample location

function p=closeEdgeLoop(p1)
  p=p1;
  step = estimateStep(p);

  % compute displacement vectors from initial sample
  dv = p;
  dv(1,:) = p(1,:)-p(1,1);
  dv(2,:) = p(2,:)-p(2,1);
  d = sqrt(sum(dv .* dv,1));

  % find the closest to the start, near the end
  kMin = closestNearEnd(d);
  if (d(kMin) > 1.7*step)
    return % not a loop
  end

  % compute dot product of step vectors with (normalized) initial step
  n = length(d);
  kLo = max([kMin-3,2]);
  kHi = min([kMin+3,n]);
  if (kHi-kLo < 1)
    return % too short for a loop
  end
  v0 = p(:,2)-p(:,1);
  v0 = v0 / norm(v0);
  for k=kLo:kHi
    d(k) = dot(v0,p(:,k)-p(:,1));
  end
  k=kLo;
  while ((d(k) < 0) && (k <= kHi))
    k=k+1;
  end
  if (k <= kLo) % not OBDVIOUSLY a loop
    return 
  end
  p(:,k) = p(:,1); % close loop
  p=p(:,1:k);
  disp(sprintf('%d/%d trimmed',n-k,n));
end

function kMin = closestNearEnd(d)
  n=length(d);
  kLo = max([n-10,2]);
  kHi = max([n-1,3]);
  dMin = d(n);
  kMin = n;
  if ((kHi-kLo<1) || (n < 3))
    return % too short an edge to look for loops
  end
  for k=kLo:kHi
    if (d(k) < dMin)
      dMin = d(k);
      kMin=k;
    end
  end
end

function step = estimateStep(p)
  n=size(p,2);
  if (n < 2)
    step = .001;
    return
  end
  dv = p(:,2:n)-p(:,1:n-1);
  d = sqrt( sum(dv .* dv,1));
  step = median(d);
end
