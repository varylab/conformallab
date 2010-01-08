package de.jtem.halfedgetools.functional.conformal.hds;

import javax.vecmath.Point3d;

import de.varylab.discreteconformal.functional.node.ConformalEdge;


public class MyConformalEdge extends ConformalEdge<MyConformalVertex, MyConformalEdge, MyConformalFace> {

   public double getLength() {
	   Point3d s = getStartVertex().getPosition();
	   Point3d t = getTargetVertex().getPosition();
	   return s.distance(t);
   }
	   
}
