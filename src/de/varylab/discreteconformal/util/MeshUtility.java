package de.varylab.discreteconformal.util;

import static java.lang.Math.abs;
import geom3d.Basis;
import geom3d.Point;
import geom3d.Vector;

import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import de.jreality.math.Matrix;
import de.jreality.math.MatrixBuilder;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.bsp.KdTree;
import de.varylab.discreteconformal.heds.bsp.KdUtility;
import de.varylab.discreteconformal.heds.bsp.KdTree.KdPosition;


public class MeshUtility {

	
	public static double[] getMinMaxCurvatureAt(
		Point p,
		double scale,
		KdTree<CoVertex> kd
	) {
		EVD evd = getTensorInformation(p, scale, kd);
		double[] eigVal= evd.getRealEigenvalues();
		LinkedList<Double> minMax = new LinkedList<Double>();
		minMax.add(eigVal[0]);
		minMax.add(eigVal[1]);
		minMax.add(eigVal[2]);
		int index = getIndexOfMinMagnitude(eigVal);
		minMax.remove(index);
		if(minMax.get(1)<minMax.get(0))
			minMax.addFirst(minMax.removeLast());
		eigVal[0]=minMax.get(0);
		eigVal[1]=minMax.get(1);
		return eigVal;
	}
	
	public static EVD getTensorInformation(
			Point p,
			double scale,
			KdTree<CoVertex> kd
	) {
		KdPosition position = new KdPosition(p);
		Collection<CoFace> faces = KdUtility.collectFacesInRadius(kd, position, scale);
		Collection<CoEdge> edges = KdUtility.collectEdgesInRadius(kd, position, scale);
		double area=0;
		for(CoFace f :faces){
			area += f.toTriangle().computeArea();
		}
		DenseMatrix matrix = new DenseMatrix(3,3);
		DenseMatrix tmp = new DenseMatrix(3,3);
		double beta = 0;
		double edgeLength = 0;
		
		for(CoEdge e : edges){
			beta = e.getCurvature();
			edgeLength = e.getLength();
			
			Vector edge = e.getVector();
			edge.normalize();
			DenseMatrix  result = new DenseMatrix(3,3);
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					result.set(i, j, edge.get()[i]*edge.get()[j]);
				}
			}
			tmp = getEdgeCurvatureTensor(e);
			matrix.add(beta*edgeLength,tmp);
		}
		matrix.scale(1/area);
		EVD evd = null;
		try {
			evd = EVD.factorize(matrix);
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
		return evd;
	}
	
	public static int getIndexOfMinMagnitude(double[] v){
		int r = 0;
		double val = abs(v[0]);
		for (int i = 1; i < v.length; i++){
			double tmp = abs(v[i]); 
			if (tmp < val) {
				r = i;
				val = tmp;
			}
		}
		return r;
	}
	
	public static Vector[] getSortedEigenVectors(EVD evd){
		Vector[] r = new Vector[]{new Vector(),new Vector(), new Vector()};
		double[] eigVal = evd.getRealEigenvalues();
		DenseMatrix eigVec = evd.getRightEigenvectors();
		double[][] eigVecArr = new double[3][3];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				eigVecArr[i][j]= eigVec.get(j,i);
			}
		}
		//get minimal magnitude
		int i3 = getIndexOfMinMagnitude(eigVal);
		int i1 = (i3 + 1) % 3;
		int i2 = (i3 + 2) % 3;
		double k1 = eigVal[i1];
		double k2 = eigVal[i2];
		r[0] = new Vector(k1 < k2 ? eigVecArr[i1] : eigVecArr[i2]);
		r[1] = new Vector(k1 < k2 ? eigVecArr[i2] : eigVecArr[i1]);
		r[2] = new Vector(eigVecArr[i3]);
		
		return r;
	}


	public static Basis getTensor(
			Point p,
			double scale,
			KdTree<CoVertex> kd
		 ){
			KdPosition position = new KdPosition(p);
			Collection<CoFace> faces = KdUtility.collectFacesInRadius(kd, position, scale);
			Collection<CoEdge> edges = KdUtility.collectEdgesInRadius(kd, position, scale);
			double area=0;
			for(CoFace f :faces){
				area += f.toTriangle().computeArea();
			}
			DenseMatrix matrix = new DenseMatrix(3,3);
			DenseMatrix tmp = new DenseMatrix(3,3);
			double beta = 0;
			double edgeLength = 0;
			
			for(CoEdge e :edges){
				beta = e.getCurvature();
				edgeLength = e.getLength();
				tmp = getEdgeCurvatureTensor(e);
				matrix.add(beta*edgeLength,tmp);
			}
			matrix.scale(1/area);
			Vector c1 = new Vector(matrix.getData()[0],matrix.getData()[1],matrix.getData()[2]);
			Vector c2 = new Vector(matrix.getData()[3],matrix.getData()[4],matrix.getData()[5]);
			Vector c3 = new Vector(matrix.getData()[6],matrix.getData()[7],matrix.getData()[8]);
			return new Basis(c1,c2,c3);
	}
	
	public static double meanEdgeLength(CoHDS mesh) {
		double result = 0.0;
		for (CoEdge e : mesh.getEdges()) {
			Point s = e.getStartVertex().getPosition();
			Point t = e.getTargetVertex().getPosition();
			result += s.distanceTo(t);
		}
		return result / mesh.numEdges() / 2;
	}
	
	public static double absoluteCurvatureAt(
			Point p,
			double scale,
			KdTree<CoVertex> kd
	){
		KdPosition position = new KdPosition(p);
		Collection<CoEdge> edges = KdUtility.collectEdgesInRadius(kd, position, scale);
		Collection<CoFace> faces = incidentFaces(edges);
		double area=0;
		for(CoFace f :faces)
			area += f.toTriangle().computeArea();
		DenseMatrix matrix = new DenseMatrix(3,3);
		DenseMatrix tmp = new DenseMatrix(3,3);
		double beta = 0;
		double edgeLength = 0;
		
		for(CoEdge e :edges){
			beta = e.getCurvature();
			edgeLength = e.getLength();
			tmp = getEdgeCurvatureTensor(e);
			matrix.add(beta*edgeLength,tmp);
		}
		matrix.scale(1/area);
		return getColumnsLength(matrix);
	}


	private static double getColumnsLength(DenseMatrix matrix) {
		Vector c1 = new Vector(matrix.getData()[0],matrix.getData()[1],matrix.getData()[2]);
		Vector c2 = new Vector(matrix.getData()[3],matrix.getData()[4],matrix.getData()[5]);
		return (abs(c1.getLength())+abs(c2.getLength()))/2.0;
	}

	private static DenseMatrix getEdgeCurvatureTensor(CoEdge e) {
		Vector edge = e.getVector();
		edge.normalize();
		DenseMatrix  result = new DenseMatrix(3,3);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				result.set(i, j, edge.get()[i]*edge.get()[j]);
			}
		}
		return result;
	}
	
	public static Collection<CoFace> incidentFaces(Collection<CoEdge> edges) {
		HashSet<CoFace> faces = new HashSet<CoFace>(edges.size() / 3);
		for (CoEdge e : edges) {
			if (e.getLeftFace() != null)
				faces.add(e.getLeftFace());
		}
		return new LinkedList<CoFace>(faces);
	}


	public static double getCurvature(CoEdge e) {
		CoFace lf = e.getLeftFace();
		CoFace rf = e.getRightFace();
		if (lf == null || rf == null)
			return 0;
		
		return curvatureSign(e)*lf.getNormal().getAngle(rf.getNormal());
	}
	
	/*
	 * 
	 * @param e, an MEdge
	 * @return -1,0,1 the sign of the angle between the left and the right face.
	 * 			negative is concave, positive if convex
	 *          
	 */
	private static double curvatureSign(CoEdge e){
		Matrix m = MatrixBuilder.euclidean().getMatrix();
		for (int i = 0; i < 3; i++) {
			m.setEntry(i, 0, e.getVector().get(i));
			m.setEntry(i, 1, e.getLeftFace().getNormal().get(i));
			m.setEntry(i, 2, e.getRightFace().getNormal().get(i));
		}
		double det = m.getDeterminant() ;
		if(Math.abs(det)< 1E-7)
			return 0;
		else if (det<0)
			return -1;
		else 
			return 1;
	}
	
}
