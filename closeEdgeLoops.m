% Close loops on a set of edges from subCannyEdgeMEX()

function ep=closeEdgeLoops(ep1)
  ep=ep1;
  for k=1:length(ep)
	  k
    ep(k).p = closeEdgeLoop(ep(k).p);
  end
end
