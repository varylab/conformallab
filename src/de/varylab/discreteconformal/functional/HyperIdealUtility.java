package de.varylab.discreteconformal.functional;

import static de.jtem.mfc.field.Complex.fromPolar;
import static de.varylab.discreteconformal.functional.Clausen.ImLi2;
import static de.varylab.discreteconformal.math.MathUtility.arcosh;
import static de.varylab.discreteconformal.math.MathUtility.arsinh;
import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static java.lang.Math.cos;
import static java.lang.Math.cosh;
import static java.lang.Math.exp;
import static java.lang.Math.sin;
import static java.lang.Math.sinh;

import java.util.List;

import org.junit.Assert;

import de.jreality.math.Matrix;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.algorithm.triangulation.Triangulator;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class HyperIdealUtility {

	/**
	 * Calculates a third length in a right angled hyperbolic hexagon
	 * with given lengths x, y, and z
	 * @return the length of the edge opposite to z
	 */
	public static double ζ_13(double x, double y, double z) {
		double cx = cosh(x);
		double cy = cosh(y);
		double cz = cosh(z);
		double sx = sinh(x);
		double sy = sinh(y);
		double n = cx*cy + cz;
		double d = sx*sy;
		return arcosh(n/d);
	}
	
	/**
	 * Calculates an edge length in a hyperbolic pentagon with one
	 * vertex at infinity
	 */
	public static double ζ_14(double x, double y) {
		double cy = cosh(y);
		double sy = sinh(y);
		double n = exp(x) + cy;
		return arcosh(n/sy);
	}
	
	/**
	 * Calculates a length in a hyperbolic quadrilateral with two
	 * vertices at infinity
	 */
	public static double ζ_15(double x) {
		return 2 * arsinh(exp(x/2));
	}
	
	/**
	 * Calculates an angle in a hyperbolic triangle with edge
	 * length x, y, and z. 
	 * @return the angle opposite of edge z.
	 */
	public static double ζ(double x, double y, double z) {
		double cx = cosh(x);
		double cy = cosh(y);
		double cz = cosh(z);
		double sx = sinh(x);
		double sy = sinh(y);
		double n = cx*cy - cz;
		double d = sx*sy;
		return acos(n/d);
	}
	
	
	/**
	 * Calculates the hyperbolic volume of the generalized 
	 * hyperbolic tetrahedron with dihedral angles α and β.
	 * Formula taken from </a href="http://dx.doi.org/10.1007/0-387-29555-0_13">here</a>.
	 * @return the volume of the tetrahedron
	 */
	public static double calculateTetrahedronVolume(double A, double B, double C, double D, double E, double F) {
		double sA = sin(A), sB = sin(B), sC = sin(C), sD = sin(D), sE = sin(E), sF = sin(F);
		double cA = cos(A), cB = cos(B), cC = cos(C), cD = cos(D), cE = cos(E), cF = cos(F);	
		Complex ad = fromPolar(1, A + D);
		Complex be = fromPolar(1, B + E);
		Complex cf = fromPolar(1, C + F);
		Complex abc = fromPolar(1, A + B + C);
		Complex abf = fromPolar(1, A + B + F);
		Complex ace = fromPolar(1, A + C + E);
		Complex aef = fromPolar(1, A + E + F);
		Complex bcd = fromPolar(1, B + C + D);
		Complex bdf = fromPolar(1, B + D + F);
		Complex def = fromPolar(1, D + E + F);
		Complex cde = fromPolar(1, C + D + E);
		Complex abde = ad.times(be);
		Complex acdf = ad.times(cf);
		Complex bcef = be.times(cf);
		Complex abcdef = abc.times(def);
		Complex z = ad.plus(be).plus(cf).plus(abf).plus(ace).plus(bcd).plus(def).plus(abcdef);
		Matrix G = new Matrix(
			1.0,	-cA,	-cB,	-cF,
			-cA,	1.0,	-cC,	-cE,
			-cB,	-cC,	1.0,	-cD,
			-cF,	-cE,	-cD,	1.0
		);
		Complex sqrtG = new Complex(G.getDeterminant()).sqrt();
		
		Complex f = new Complex(sA*sD + sB*sE + sC*sF);
		Complex f1 = f.minus(sqrtG);
		Complex f2 = f.plus(sqrtG);
		Complex z1 = f1.times(-2).divide(z);
		Complex z2 = f2.times(-2).divide(z);
		
		double U1 = 0.5 * (
			+ ImLi2(z1) 
			+ ImLi2(abde.times(z1)) 
			+ ImLi2(acdf.times(z1)) 
			+ ImLi2(bcef.times(z1)) 
			- ImLi2(abc.times(z1).times(-1)) 
			- ImLi2(aef.times(z1).times(-1)) 
			- ImLi2(bdf.times(z1).times(-1)) 
			- ImLi2(cde.times(z1).times(-1))
		); 
		double U2 = 0.5 * (
			+ ImLi2(z2) 
			+ ImLi2(abde.times(z2)) 
			+ ImLi2(acdf.times(z2)) 
			+ ImLi2(bcef.times(z2)) 
			- ImLi2(abc.times(z2).times(-1)) 
			- ImLi2(aef.times(z2).times(-1)) 
			- ImLi2(bdf.times(z2).times(-1)) 
			- ImLi2(cde.times(z2).times(-1))
		);
		
		return (U1 - U2) / 2;
	}
	
	
	public static CoHDS createLawsonsSurfaceWithBranchPoints() {
		CoHDS hds = createLawsonsSurface();
//		StellarLinear sl = new StellarLinear();
		
		return hds;
	}
	

	public static CoHDS createLawsonsSurface() {
		CoHDS hds = new CoHDS(); 
		CoVertex A = hds.addNewVertex();
		CoVertex B = hds.addNewVertex();
		CoVertex C = hds.addNewVertex();
		CoVertex D = hds.addNewVertex();
		hds.addNewEdges(24);
		hds.getEdge(0).linkOppositeEdge(hds.getEdge(6));
		hds.getEdge(1).linkOppositeEdge(hds.getEdge(21));
		hds.getEdge(2).linkOppositeEdge(hds.getEdge(4));
		hds.getEdge(3).linkOppositeEdge(hds.getEdge(23));
		hds.getEdge(5).linkOppositeEdge(hds.getEdge(11));
		hds.getEdge(7).linkOppositeEdge(hds.getEdge(9));
		hds.getEdge(8).linkOppositeEdge(hds.getEdge(14));
		hds.getEdge(10).linkOppositeEdge(hds.getEdge(12));
		hds.getEdge(13).linkOppositeEdge(hds.getEdge(19));
		hds.getEdge(15).linkOppositeEdge(hds.getEdge(17));
		hds.getEdge(16).linkOppositeEdge(hds.getEdge(22));
		hds.getEdge(18).linkOppositeEdge(hds.getEdge(20));
		hds.getEdge(0).setTargetVertex(A);
		hds.getEdge(1).setTargetVertex(B);
		hds.getEdge(2).setTargetVertex(D);
		hds.getEdge(3).setTargetVertex(C);
		hds.getEdge(4).setTargetVertex(B);
		hds.getEdge(5).setTargetVertex(A);
		hds.getEdge(6).setTargetVertex(C);
		hds.getEdge(7).setTargetVertex(D);
		hds.getEdge(8).setTargetVertex(D);
		hds.getEdge(9).setTargetVertex(C);
		hds.getEdge(10).setTargetVertex(A);
		hds.getEdge(11).setTargetVertex(B);
		hds.getEdge(12).setTargetVertex(C);
		hds.getEdge(13).setTargetVertex(D);
		hds.getEdge(14).setTargetVertex(B);
		hds.getEdge(15).setTargetVertex(A);
		hds.getEdge(16).setTargetVertex(A);
		hds.getEdge(17).setTargetVertex(B);
		hds.getEdge(18).setTargetVertex(D);
		hds.getEdge(19).setTargetVertex(C);
		hds.getEdge(20).setTargetVertex(B);
		hds.getEdge(21).setTargetVertex(A);
		hds.getEdge(22).setTargetVertex(C);
		hds.getEdge(23).setTargetVertex(D);
		
		hds.getEdge(0).linkNextEdge(hds.getEdge(1));
		hds.getEdge(1).linkNextEdge(hds.getEdge(2));
		hds.getEdge(2).linkNextEdge(hds.getEdge(3));
		hds.getEdge(3).linkNextEdge(hds.getEdge(0));
		
		hds.getEdge(4).linkNextEdge(hds.getEdge(5));
		hds.getEdge(5).linkNextEdge(hds.getEdge(6));
		hds.getEdge(6).linkNextEdge(hds.getEdge(7));
		hds.getEdge(7).linkNextEdge(hds.getEdge(4));
		
		hds.getEdge(8).linkNextEdge(hds.getEdge(9));
		hds.getEdge(9).linkNextEdge(hds.getEdge(10));
		hds.getEdge(10).linkNextEdge(hds.getEdge(11));
		hds.getEdge(11).linkNextEdge(hds.getEdge(8));
		
		hds.getEdge(12).linkNextEdge(hds.getEdge(13));
		hds.getEdge(13).linkNextEdge(hds.getEdge(14));
		hds.getEdge(14).linkNextEdge(hds.getEdge(15));
		hds.getEdge(15).linkNextEdge(hds.getEdge(12));
		
		hds.getEdge(16).linkNextEdge(hds.getEdge(17));
		hds.getEdge(17).linkNextEdge(hds.getEdge(18));
		hds.getEdge(18).linkNextEdge(hds.getEdge(19));
		hds.getEdge(19).linkNextEdge(hds.getEdge(16));
		
		hds.getEdge(20).linkNextEdge(hds.getEdge(21));
		hds.getEdge(21).linkNextEdge(hds.getEdge(22));
		hds.getEdge(22).linkNextEdge(hds.getEdge(23));
		hds.getEdge(23).linkNextEdge(hds.getEdge(20));
		
		HalfEdgeUtils.fillAllHoles(hds);
		List<CoEdge> auxEdges = Triangulator.triangulateSingleSource(hds);
		
		Assert.assertTrue(HalfEdgeUtils.isValidSurface(hds, true));
		Assert.assertEquals(2, HalfEdgeUtils.getGenus(hds));
		
		
		int index = 0;
		for (CoVertex v : hds.getVertices()) {
			v.setSolverIndex(index++);
			v.setTheta(2 * PI);
		}
		for (CoEdge e : hds.getPositiveEdges()) {
			e.setSolverIndex(index);
			e.getOppositeEdge().setSolverIndex(index++);
			if (auxEdges.contains(e)) {
				e.setTheta(PI);
				e.getOppositeEdge().setTheta(PI);
			} else {
				e.setTheta(PI/2);
				e.getOppositeEdge().setTheta(PI/2);
			}
		}
		return hds;
	}
	
}
