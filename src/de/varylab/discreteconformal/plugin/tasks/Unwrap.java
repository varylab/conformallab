package de.varylab.discreteconformal.plugin.tasks;

import static de.varylab.discreteconformal.util.CuttingUtility.cutManifoldToDisk;

import javax.swing.SwingWorker;

import no.uib.cipr.matrix.Vector;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.plugin.adapter.HyperbolicLengthWeightAdapter;
import de.varylab.discreteconformal.unwrapper.CHyperbolicLayout;
import de.varylab.discreteconformal.unwrapper.CHyperbolicUnwrapper;
import de.varylab.discreteconformal.unwrapper.CHyperbolicUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.Genus0Layout;
import de.varylab.discreteconformal.unwrapper.Genus0Unwrapper;
import de.varylab.discreteconformal.unwrapper.Genus0UnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.Unwrapper;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class Unwrap extends SwingWorker<CoHDS, Object> {

	private CoHDS
		surface = null;
	private boolean
		usePetsc = false,
		quantizeCones = false;
	private int
		numCones = 0;
	
	// result values -------
	public int 
		genus = -1;
	public CoVertex
		layoutRoot = null;
	public CuttingInfo<CoVertex, CoEdge, CoFace> 
		cutInfo = null;
	
	
	public Unwrap(CoHDS surface) {
		this.surface = surface;
	}
	
	
	@Override
	protected CoHDS doInBackground() throws Exception {
		setProgress(0);
		genus = HalfEdgeUtils.getGenus(surface);
		Vector u = null;
		switch (genus) {
		// disk or torus ---------------------
		case 0: 
			Unwrapper unwrapper = null;
			if (usePetsc) {
				unwrapper = new Genus0UnwrapperPETSc();
			} else {
				unwrapper = new Genus0Unwrapper();
			}
			u = unwrapper.unwrap(surface, numCones, quantizeCones);
			setProgress(50);
			layoutRoot = Genus0Layout.doLayout(surface, u);
			setProgress(100);
			return surface;
		// sphere ----------------------------
		case 1:
			return null;
		// genus > 1 -------------------------
		default:
			if (usePetsc) {
				unwrapper = new CHyperbolicUnwrapperPETSc();
			} else {
				unwrapper = new CHyperbolicUnwrapper();
			}
			u = unwrapper.unwrap(surface, 0, false);
			setProgress(50);
			HyperbolicLengthWeightAdapter hypWa = new HyperbolicLengthWeightAdapter(u);
			CoVertex root = surface.getVertex(getMinUIndex(u));
			cutInfo = cutManifoldToDisk(surface, root, hypWa);
			layoutRoot = CHyperbolicLayout.doLayout(surface, u);
			setProgress(100);
			return surface;
		}
	}
	
	
	private int getMinUIndex(Vector u) {
		int index = 0;
		double iVal = u.get(0);
		for (int i = 1; i < u.size(); i++) {
			double val = u.get(i);
			if (iVal < val) {
				index = i;
				iVal = i;
			}
		}
		return index;
	}


	public CoHDS getSurface() {
		return surface;
	}


	public void setSurface(CoHDS surface) {
		this.surface = surface;
	}


	public boolean isUsePetsc() {
		return usePetsc;
	}


	public void setUsePetsc(boolean usePetsc) {
		this.usePetsc = usePetsc;
	}


	public boolean isQuantizeCones() {
		return quantizeCones;
	}


	public void setQuantizeCones(boolean quantizeCones) {
		this.quantizeCones = quantizeCones;
	}


	public int getNumCones() {
		return numCones;
	}


	public void setNumCones(int numCones) {
		this.numCones = numCones;
	}
	

}
