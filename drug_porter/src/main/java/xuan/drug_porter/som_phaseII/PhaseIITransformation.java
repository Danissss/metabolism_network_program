package xuan.drug_porter.som_phaseII;

import java.util.List;

import org.openscience.cdk.Atom;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import xuan.drug_porter.som_phaseII.FunctionalGroup;

public class PhaseIITransformation {

static IAtom oxygen = new Atom(8);
	
	public static FunctionalGroup FG = new FunctionalGroup();
	
	/**
	 * do the transformation by following the rule system
	 * @param mole
	 * @param som
	 * @return
	 */
	public static IAtomContainer phaseItransformer(IAtomContainer mole, int som) {
		
		int num_of_atom = mole.getAtomCount();
		int num_of_bond = mole.getBondCount();
		
		IAtom site_of_atom = mole.getAtom(som);
		List<IAtom> connected_atom = mole.getConnectedAtomsList(site_of_atom);
		for(int i = 0; i < connected_atom.size(); i++) {
			
			// rule system here: add_enzyme;
			
		}
		
		
		
		
		
		return mole;
	}
	
	
	
	
	public static void main(String[] args) {
		
		
		
		
	}
}
