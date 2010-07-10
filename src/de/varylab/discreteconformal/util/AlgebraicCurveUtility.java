package de.varylab.discreteconformal.util;

import java.util.Set;

import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class AlgebraicCurveUtility {

	
	public static void calculateModulus(
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo, 
		CoVertex root
	) {
		Set<CoVertex> copies = cutInfo.getCopies(root);
		System.out.println(copies);
	}
	
	
}
