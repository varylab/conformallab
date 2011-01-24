package de.varylab.discreteconformal.uniformization;

import static java.math.BigDecimal.ONE;
import static java.math.BigDecimal.ZERO;

import java.math.BigDecimal;

import de.jreality.math.Matrix;

public class FundamentalEdge implements Comparable<FundamentalEdge> {
	
	public int 
		// index in the clockwise ordering
		index = 0, 
		// number of edges in the cut mesh
		sourceEdgeCount = 0;
	public FundamentalVertex 
		start = new FundamentalVertex(0),
		end = new FundamentalVertex(0);
	public Matrix
		// identification motion
		motion = new Matrix(); 
	public BigDecimal[]
	    startPosition = {ZERO, ZERO, ZERO, ONE},
	    motionBig = null;
	public FundamentalEdge
		prevEdge = null,
		nextEdge = null;
	public FundamentalEdge
		partner = null;

	public FundamentalEdge(int index) {
		this.index = index;
	}

	public FundamentalEdge(
		int index,
		int sourceEdges,
		FundamentalVertex start, 
		FundamentalVertex end,
		BigDecimal[] startPosition,
		Matrix motion,
		BigDecimal[] motionBig
	) {
		this(index);
		this.sourceEdgeCount = sourceEdges;
		this.start = start;
		this.end = end;
		this.startPosition = startPosition;
		this.motion = motion;
		this.motionBig = motionBig;
	}
	
	@Override
	public String toString() {
		return "FunEdge " + index;
	}

	@Override
	public int compareTo(FundamentalEdge o) {
		return index - o.index;
	}
	
}