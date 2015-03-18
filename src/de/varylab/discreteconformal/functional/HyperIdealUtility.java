package de.varylab.discreteconformal.functional;

import static de.jtem.mfc.field.Complex.fromPolar;
import static de.varylab.discreteconformal.functional.Clausen.Л;
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
import de.jreality.math.Matrix;
import de.jtem.mfc.field.Complex;


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
	 * Calculates the hyperbolic volume of the hyperidel tetrahedron with dihedral angles
	 * γ1, γ2, and γ3 at the ideal vertex. Formula taken from <a href="http://arxiv.org/abs/math/0603097">here</a>
	 * @return the volume of the hyperideal tetrahedron
	 */
	public static double calculateTetrahedronVolumeWithIdealVertexAtGamma(
		double γ1, double γ2, double γ3,
		double α12, double α23, double α31 
	) {
		double result = Л(γ1) + Л(γ2) + Л(γ3);
		result += Л((PI + α31 - α12 - γ1)/2) + Л((PI + α12 - α23 - γ2)/2) - Л((PI + α23 - α31 - γ3)/2);
		result += Л((PI - α31 + α12 - γ1)/2) + Л((PI - α12 + α23 - γ2)/2) - Л((PI - α23 + α31 - γ3)/2);
		result += Л((PI + α31 + α12 - γ1)/2) + Л((PI + α12 + α23 - γ2)/2) - Л((PI + α23 + α31 - γ3)/2);
		result += Л((PI - α31 - α12 - γ1)/2) + Л((PI - α12 - α23 - γ2)/2) - Л((PI - α23 - α31 - γ3)/2);
		return result / 2;
	}
	
	/**
	 * Calculates the hyperbolic volume of the generalized 
	 * hyperbolic tetrahedron with dihedral angles α and β.
	 * Formula taken from <a href="http://dx.doi.org/10.1007/0-387-29555-0_13">here</a>.
	 * @return the volume of the tetrahedron
	 */
	public static double calculateTetrahedronVolume(double A, double B, double C, double D, double E, double F) {
		int piCounter = 0;
		piCounter += A == PI ? 1 : 0;
		piCounter += B == PI ? 1 : 0;
		piCounter += C == PI ? 1 : 0;
		piCounter += D == PI ? 1 : 0;
		piCounter += E == PI ? 1 : 0;
		piCounter += F == PI ? 1 : 0;
		if (piCounter > 0) return 0.0;  // degenerate tetrahedron
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
	
}
