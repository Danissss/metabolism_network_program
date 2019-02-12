package xuan.drug_porter.som_phaseI;

import java.util.List;

import org.openscience.cdk.Atom;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import xuan.drug_porter.descriptorUtils.SdfToSample;
import xuan.drug_porter.som_phaseI.SomPrediction;



public class PhaseITransformation {
	
	
	
	static IAtom oxygen = new Atom(8);
	
	
	
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
			
			// rule system here:
			
		}
		
		
		
		
		
		return mole;
	}
	
	
	/**
	 * this is for hydroxylation reaction
	 * @return
	 */
	public static IAtomContainer add_OH_group(IAtomContainer mole, int som) {
		
		IAtomContainer enzyme_mole_no_h = AtomContainerManipulator.removeHydrogens(mole);
		int num_of_atoms = enzyme_mole_no_h.getAtomCount();
		int new_atom_index = num_of_atoms - 1 + 1;
		enzyme_mole_no_h.addAtom(oxygen);
		enzyme_mole_no_h.addBond(som, new_atom_index, IBond.Order.SINGLE);
	
		return enzyme_mole_no_h;
		
	}
	
	
	
	public static void main(String[] args) {
		
		
		
		
	}

}
