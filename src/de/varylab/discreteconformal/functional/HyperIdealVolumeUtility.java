package de.varylab.discreteconformal.functional;

import static de.jtem.mfc.field.Complex.fromPolar;
import static de.varylab.discreteconformal.functional.Clausen.ImLi2;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import de.jreality.math.Matrix;
import de.jtem.mfc.field.Complex;

public class HyperIdealVolumeUtility {

	
	/**
	 * Calculate the hyperbolic volume of the generalized 
	 * hyperbolic tetrahedron with dihedral angles α and β.
	 * Taken from 
	 * @return the volume of the tetrahedron
	 */
	public static double calculateVolume(double A, double B, double C, double D, double E, double F) {
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