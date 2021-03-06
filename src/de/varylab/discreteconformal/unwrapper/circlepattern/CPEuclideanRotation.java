/**
This file is part of a jTEM project.
All jTEM projects are licensed under the FreeBSD license 
or 2-clause BSD license (see http://www.opensource.org/licenses/bsd-license.php). 

Copyright (c) 2002-2010, Technische Universität Berlin, jTEM
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

-	Redistributions of source code must retain the above copyright notice, 
	this list of conditions and the following disclaimer.

-	Redistributions in binary form must reproduce the above copyright notice, 
	this list of conditions and the following disclaimer in the documentation 
	and/or other materials provided with the distribution.
 
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS 
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY 
OF SUCH DAMAGE.
**/

package de.varylab.discreteconformal.unwrapper.circlepattern;

import static java.lang.Math.PI;

import java.util.Map;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Vertex;
import de.jtem.mfc.field.Complex;
import de.jtem.mfc.geometry.ComplexProjective1;
import de.jtem.mfc.group.Moebius;
import de.varylab.discreteconformal.unwrapper.circlepattern.CirclePatternLayout.Rotation;

public class CPEuclideanRotation 
<
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>
> implements Rotation<V, E, F> {

	
	@Override
	public double[] rotate(double[] p, double[] center, Double phi, Double logScale) {
		Moebius rot = new Moebius();
		Complex c1 = new Complex(center[0], center[1]);
		Complex c2 = new Complex(1,0);
		ComplexProjective1 c = new ComplexProjective1(c1, c2);
		rot.assignEuclideanLogScaleRotation(c, logScale, phi);
		Complex pc = new Complex(p[0], p[1]);
		pc = rot.applyTo(pc);
		return new double[] {pc.re, pc.im};
	}

	
	@Override
	public double getPhi(E edge, Map<F, Double> rhoMap, Map<E, Double> thetaMap) {
		double th = thetaMap.get(edge);
		double thStar = PI - th;
		if (edge.getLeftFace() == null || edge.getRightFace() == null)
			return thStar;
		double leftRho = rhoMap.get(edge.getLeftFace());
		double rightRho = rhoMap.get(edge.getRightFace());
		double p = p(thStar, rightRho - leftRho);
		return 0.5*(p + thStar);
	}

	
	@Override
	public Double getRadius(Double rho) {
		return Math.exp(rho);
	}
	
    public Double p(Double thStar, Double diffRho) {
        Double exp = Math.exp(diffRho);
        Double tanhDiffRhoHalf = (exp - 1.0) / (exp + 1.0);
        return  2.0 * Math.atan(Math.tan(0.5 * thStar) * tanhDiffRhoHalf);
    }
	
	
}
