package de.varylab.discreteconformal;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;


public class ConformalAdapterSet extends AdapterSet {

	private static final long 
		serialVersionUID = 1L;

	public ConformalAdapterSet() {
		add(new CoPositionAdapter());
		add(new CoTexturePositionAdapter());
		addAll(AdapterSet.createGenericAdapters());
	}
	
}  
