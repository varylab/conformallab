package de.varylab.discreteconformal.unwrapper;

import java.util.Map;

import de.jtem.halfedgetools.adapter.AdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public interface Unwrapper {

	public void unwrap(CoHDS surface, int g, AdapterSet aSet) throws Exception;
	
	public void setGradientTolerance(double tol);
	public void setMaxIterations(int maxIterations);
	public void setCutRoot(CoVertex root);
	
	public CuttingInfo<CoVertex, CoEdge, CoFace> getCutInfo();
	public Map<CoEdge, Double> getlengthMap();
	public CoVertex getLayoutRoot();
	
}
