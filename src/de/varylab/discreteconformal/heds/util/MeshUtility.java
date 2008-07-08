package de.varylab.discreteconformal.heds.util;

import static java.lang.Math.abs;
import static java.lang.Math.min;
import geom3d.Basis;
import geom3d.Point;
import geom3d.Quad;
import geom3d.Vector;

import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import de.varylab.discreteconformal.heds.CEdge;
import de.varylab.discreteconformal.heds.CFace;
import de.varylab.discreteconformal.heds.CVertex;
import de.varylab.discreteconformal.heds.HDS;
import de.varylab.discreteconformal.heds.bsp.KdTree;
import de.varylab.discreteconformal.heds.bsp.KdUtility;
import de.varylab.discreteconformal.heds.bsp.KdTree.KdPosition;


public class MeshUtility {

	
	public static double[] getMinMaxCurvatureAt(
		Point p,
		double scale,
		KdTree<CVertex> kd
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
			KdTree<CVertex> kd
	) {
		KdPosition position = new KdPosition(p);
		Collection<CFace> faces = KdUtility.collectFacesInRadius(kd, position, scale);
		Collection<CEdge> edges = KdUtility.collectEdgesInRadius(kd, position, scale);
		double area=0;
		for(CFace f :faces){
			area += f.toTriangle().computeArea();
		}
		DenseMatrix matrix = new DenseMatrix(3,3);
		DenseMatrix tmp = new DenseMatrix(3,3);
		double beta = 0;
		double edgeLength = 0;
		
		for(CEdge e : edges){
			beta = e.getAngle();
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
			KdTree<CVertex> kd
		 ){
			KdPosition position = new KdPosition(p);
			Collection<CFace> faces = KdUtility.collectFacesInRadius(kd, position, scale);
			Collection<CEdge> edges = KdUtility.collectEdgesInRadius(kd, position, scale);
			double area=0;
			for(CFace f :faces){
				area += f.toTriangle().computeArea();
			}
			DenseMatrix matrix = new DenseMatrix(3,3);
			DenseMatrix tmp = new DenseMatrix(3,3);
			double beta = 0;
			double edgeLength = 0;
			
			for(CEdge e :edges){
				beta = e.getAngle();
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
	
	public static double meanEdgeLength(HDS mesh) {
		double result = 0.0;
		for (CEdge e : mesh.getEdges()) {
			Point s = e.getStartVertex().getPosition();
			Point t = e.getTargetVertex().getPosition();
			result += s.distanceTo(t);
		}
		return result / mesh.numEdges() / 2;
	}
	
	public static double absoluteCurvatureAt(
			Point p,
			double scale,
			KdTree<CVertex> kd
	){
		KdPosition position = new KdPosition(p);
		Collection<CEdge> edges = KdUtility.collectEdgesInRadius(kd, position, scale);
		Collection<CFace> faces = incidentFaces(edges);
		double area=0;
		for(CFace f :faces)
			area += f.toTriangle().computeArea();
		DenseMatrix matrix = new DenseMatrix(3,3);
		DenseMatrix tmp = new DenseMatrix(3,3);
		double beta = 0;
		double edgeLength = 0;
		
		for(CEdge e :edges){
			beta = e.getAngle();
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

	private static DenseMatrix getEdgeCurvatureTensor(CEdge e) {
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
	
	public static Collection<CFace> incidentFaces(Collection<CEdge> edges) {
		HashSet<CFace> faces = new HashSet<CFace>(edges.size() / 3);
		for (CEdge e : edges) {
			if (e.getLeftFace() != null)
				faces.add(e.getLeftFace());
		}
		return new LinkedList<CFace>(faces);
	}
	
	public static double getSpacingDistance(double k, double eps){
		return 2*Math.sqrt(eps*(2/abs(k)-eps));
	}
	
	public static double getSpacingDistance(double k){
		double eps = 0.008;
		return 2*Math.sqrt(eps*(2/abs(k)-eps));
	}
	
	public static Quad createQuad(Basis basis, Point center,double kMin,double kMax, double scale ){
		Quad quad = new Quad();
		double distMin = getSpacingDistance(kMin);
		Vector dMin = new Vector(basis.getX()).scaleTo(distMin);
		
		double distMax = getSpacingDistance(kMax);
		Vector dMax = new Vector(basis.getY()).scaleTo(distMax);
		if(distMin>scale*20)
			dMin = new Vector(basis.getX()).scaleTo(scale*50);
		
		if (distMax>scale*20)
			dMax = new Vector(basis.getY()).scaleTo(scale*50);
		
		Point min1 = new Point(center).add(dMin).asPoint();
		
		quad.setA( new Point(min1).add(dMax).asPoint());
		quad.setB( new Point(min1).add(dMax.times(-1)).asPoint());
		Point min2 = new Point(center).add(dMin.times(-1)).asPoint();
		quad.setC( new Point(min2).add(dMax).asPoint());
		quad.setD(new Point(min2).add(dMax.times(-1)).asPoint());

		return quad;
	}
	public static Quad createQuad(Basis basis, Point center,double kMin,double kMax,double scale, double maxScale ){
		Quad quad = new Quad();
		
		double distMin = getSpacingDistance(kMin, scale);
		distMin = min(distMin, maxScale);
		Vector dMin = new Vector(basis.getX()).scaleTo(distMin);
		
		double distMax = getSpacingDistance(kMax, scale);
		distMax = min(distMax, maxScale);
		Vector dMax = new Vector(basis.getY()).scaleTo(distMax);
		
		Point min1 = new Point(center).add(dMin).asPoint();
		
		quad.setA( new Point(min1).add(dMax).asPoint());
		quad.setB( new Point(min1).add(dMax.times(-1)).asPoint());
		Point min2 = new Point(center).add(dMin.times(-1)).asPoint();
		quad.setC( new Point(min2).add(dMax).asPoint());
		quad.setD(new Point(min2).add(dMax.times(-1)).asPoint());

		return quad;
	}
//	test usage
//	public static void main(String[] args) {
//		Basis basis = new Basis();
//		basis.setX(new Vector(1,0,0));
//		basis.setY(new Vector(0,1,0));
//		basis.setZ(new Vector(0,0,1));
//		Point center = new Point(0,0,0);
//		Quad result = createQuad(basis, center, 0.01, 0.1);
//		System.err.println(result.toString());
//	}
}
