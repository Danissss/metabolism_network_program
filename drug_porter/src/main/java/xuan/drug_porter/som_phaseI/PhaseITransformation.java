package xuan.drug_porter.som_phaseI;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.math.plot.utils.Array;
import org.openscience.cdk.Atom;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.modeling.builder3d.ModelBuilder3D;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import weka.core.Instances;
import xuan.drug_porter.descriptorUtils.SdfToSample;
import xuan.drug_porter.som_phaseI.SomPrediction;



public class PhaseITransformation {
	
	
	
	static IAtom oxygen = new Atom(8);
	public static SmilesGenerator smigen = new SmilesGenerator(SmiFlavor.Isomeric);
	
	
	public static IAtomContainer Get3DAtomContainerFromSmiles(String smiles) throws CDKException, CloneNotSupportedException, IOException {
		
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer original_mole = builder.newInstance(IAtomContainer.class);
		IAtomContainer molecule_3D = builder.newInstance(IAtomContainer.class);
		SmilesParser temp_smiles = new SmilesParser(builder);
	 	IAtomContainer atom_container = temp_smiles.parseSmiles(smiles);
	 	AtomContainerManipulator.suppressHydrogens(atom_container);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(atom_container);
	
	 	StructureDiagramGenerator sdg = new StructureDiagramGenerator();
		sdg.setMolecule(atom_container);
		sdg.generateCoordinates();
		original_mole = sdg.getMolecule();
		
		
		ModelBuilder3D mb3d = ModelBuilder3D.getInstance(builder);
		molecule_3D = mb3d.generate3DCoordinates(original_mole, false);
		return molecule_3D;
	}
	
	/**
	 * do the transformation by following the rule system
	 * may return two metabolites due to cleavge; may return one metabolites by hydroxylation
	 * @param mole
	 * @param som
	 * @return
	 * @throws CDKException 
	 * @throws IOException 
	 * @throws CloneNotSupportedException 
	 */
	public static ArrayList<IAtomContainer> phaseItransformer(IAtomContainer mole, int som) throws CDKException, CloneNotSupportedException, IOException {
		
		
		ArrayList<IAtomContainer> return_mole = new ArrayList<IAtomContainer>(); 
		int num_of_atom = mole.getAtomCount();
		int num_of_bond = mole.getBondCount();
		
		IAtom site_of_atom = mole.getAtom(som);
		List<IAtom> connected_atom = mole.getConnectedAtomsList(site_of_atom);
		List<IBond>	connected_bond = mole.getConnectedBondsList(site_of_atom);
		

		
		
		ArrayList<String> symbol_pool = new ArrayList<String>();
		int connected_atom_size = connected_atom.size();
		for(int i = 0; i < connected_atom_size; i++) {
			
			IAtom connected_atom_obj = connected_atom.get(i);
			symbol_pool.add(connected_atom_obj.getSymbol());
//			System.out.println(Integer.toString(connected_atom_obj.getIndex()) + ":" +connected_atom_obj.getSymbol());
			
		}
		
		ArrayList<IBond.Order> bond_pool = new ArrayList<IBond.Order>();
		for(int i = 0; i < connected_bond.size(); i++) {
			IBond bond = connected_bond.get(i);
			bond_pool.add(bond.getOrder());
		}
		
		
		// SOM is Carbon
		if(site_of_atom.getSymbol() == "C") {
			// if the site_of_atom is carbon
			// hydroxylation
			// epoxidation
			// note: all the cleavge the site_of_metabolism is on Carbon
			// C-Oxidation: C-Oxidation might comes from Hydroxylation, because of unstable electron field of oxygen, single bond becomes double bond
			//
			
			if(connected_atom.size() == 1) {
				
				if(symbol_pool.contains("O")) {
//					System.out.println("O-demethylation");
					// O-demethylation
					IAtom other_atom = connected_atom.get(0);
					return_mole = cleavage(mole,site_of_atom,other_atom);
					
					
				}else if(symbol_pool.contains("N")) {
					// N-demethylation
					IAtom other_atom = connected_atom.get(0);
					return_mole = cleavage(mole,site_of_atom,other_atom);
					
					
				}else if(symbol_pool.contains("C")) {
					// aliphatic hydroxylation
					IAtomContainer biotranformered = add_OH_group(mole,som);
					return_mole.add(biotranformered);
				}
				
				
			}
			else if(connected_atom.size() == 2) {
				// one or more methyl group
				if(symbol_pool.contains("C") && symbol_pool.contains("O")) {
					// O-demethylation with potentially more than one methyl group
					
					for(int i = 0; i< connected_atom_size; i++) {
						if(connected_atom.get(i).getSymbol() == "O") {
							IAtom other_atom = connected_atom.get(i);
							return_mole = cleavage(mole,site_of_atom,other_atom);
						}
					}
					
					
				}else if(symbol_pool.contains("C") && symbol_pool.contains("N")) {
					// N-demethylation with potentially more than one methyl group
					// some wired case may happen (for example: CCN) but let's hope that our machine learning model doesn't be that stupid...
					for(int i = 0; i< connected_atom_size; i++) {
						if(connected_atom.get(i).getSymbol() == "N") {
							IAtom other_atom = connected_atom.get(i);
							return_mole = cleavage(mole,site_of_atom,other_atom);
						}
					}
					
					
				}else if(symbol_pool.contains("C") && !symbol_pool.contains("N") && !symbol_pool.contains("O") && !symbol_pool.contains("S")) {
					// check if it is aromatic hydroxylation or aliphatic hydroxylation or epoxidation
					IAtomContainer biotranformered = add_OH_group(mole,som);
					return_mole.add(biotranformered);
					// epoxidation may be hard to capture
					
					
					
				}
				
				
			}else if(connected_atom.size() == 3) {
				
				
			}
			
			
			
		}
		// SOM is Oxygen
		else if(site_of_atom.getSymbol() == "O") {
			// if the site_of_atom is oxygen; then two possible reaction: Reduction
			// Desulfation
			if(symbol_pool.contains("C") && bond_pool.contains(IBond.Order.DOUBLE)) {
				// reduction
				IAtomContainer biotranformered = OxygenReduction(mole,som);
			}
			
		}
		
		// SOM is Sulfur
		else if(site_of_atom.getSymbol() == "S") {
			// if predictor point to S atom, it probably is S-Oxidation
			// s-reduction
//			if(connected_atom) {
//				
//			}else
//			IAtomContainer biotransformed = Oxidation(mole,som);
//			return_mole.add(biotransformed);
			
		}
		 
		// SOM is Nitrogen
		else if(site_of_atom.getSymbol() == "N") {
			IAtomContainer biotransformed = Oxidation(mole,som);
			return_mole.add(biotransformed);
			
		}
		
		
		

		
		
		
		
		// reduction on C=O -> C-OH; NO2 -> NH2; 
		// reduction som is on the oxygen
		
		return return_mole;
	}
	
	
	/**
	 * perform cleavage of bond of two atom;
	 * @param mole
	 * @param atom1
	 * @param atom2
	 * @return
	 * @throws CDKException 
	 * @throws IOException 
	 * @throws CloneNotSupportedException 
	 */
	public static ArrayList<IAtomContainer> cleavage(IAtomContainer mole, IAtom atom1, IAtom atom2) throws CDKException, CloneNotSupportedException, IOException {
		int ignore_metabolites_factor = 4;
		ArrayList<IAtomContainer> atom_list = new ArrayList<IAtomContainer>();
		
		IBond bonds = mole.getBond(atom1, atom2);
		mole.removeBond(bonds);
		String smiles = smigen.create(mole);
		if(smiles.contains(".")) {
			String[] smiles_list = smiles.split("\\.");
			for(int i = 0; i < smiles_list.length; i++) {
				IAtomContainer tmp_container = Get3DAtomContainerFromSmiles(smiles_list[i]);
				IAtomContainer tmp_container_no_h = AtomContainerManipulator.removeHydrogens(tmp_container);
				if(tmp_container_no_h.getAtomCount() > ignore_metabolites_factor) {
					atom_list.add(tmp_container_no_h);
				}
			}
		}
		
		
		return atom_list;
	}
	
	
	
	/**
	 * (require test)
	 * reduction for Oxygen
	 * @param mole
	 * @param som
	 * @return
	 */
	public static IAtomContainer OxygenReduction(IAtomContainer mole, int som) {
		
		IAtomContainer enzyme_mole_no_h = AtomContainerManipulator.removeHydrogens(mole);
		IAtom site_of_atom = enzyme_mole_no_h.getAtom(som);
//		List<IAtom> connected_atom = enzyme_mole_no_h.getConnectedAtomsList(site_of_atom);
		List<IBond>	connected_bond = enzyme_mole_no_h.getConnectedBondsList(site_of_atom);
		
		for(int i = 0; i < connected_bond.size(); i++) {
			IBond bond = connected_bond.get(i);
			IAtom other_atoms = bond.getOther(site_of_atom);
			if(bond.getOrder() == IBond.Order.DOUBLE && other_atoms.getSymbol() == "C") {
				// reduction
				bond.setOrder(IBond.Order.SINGLE);
			}
		}
		
		return enzyme_mole_no_h;
		
		
	}
	
	/**
	 * oxidation for S-Oxidation and N-Oxidation;
	 * @param mole
	 * @param som
	 * @return
	 * @throws CDKException
	 * @throws CloneNotSupportedException
	 * @throws IOException
	 */
	public static IAtomContainer Oxidation(IAtomContainer mole, int som) throws CDKException, CloneNotSupportedException, IOException {
		
		
		IAtomContainer enzyme_mole_no_h = AtomContainerManipulator.removeHydrogens(mole);
		int num_of_atoms = enzyme_mole_no_h.getAtomCount();
		int new_atom_index = num_of_atoms - 1 + 1;
		enzyme_mole_no_h.addAtom(oxygen);
		enzyme_mole_no_h.addBond(som, new_atom_index, IBond.Order.DOUBLE);		
		
		return enzyme_mole_no_h;
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
	
	
	
	public static void main(String[] args) throws Exception {
		String test_smiles = "/Users/xuan/Desktop/HMDB00001.sdf";
		String smile = "C=12C(=C(OC(C([H])([H])[H])([H])[H])C(=C(C1C(C(Cl)([H])[H])=C(C(=O)O2)[H])[H])[H])[H]";   // O-demethylation
		
		IAtomContainer original_mole = Get3DAtomContainerFromSmiles(smile);
		IAtomContainer original_mole_no_h = AtomContainerManipulator.removeHydrogens(original_mole);
		
		SomPrediction sp = new SomPrediction();
		Instances ins = SomPrediction.create_test_instance(smile);
		HashMap<Integer,String> result = SomPrediction.runSomClassifier(ins, "CYP1A2");
		for (Integer key : result.keySet()) {
//			System.out.println(Integer.toString(key)+" : "+ result.get(key));
			if(result.get(key) == "Yes") {
				ArrayList<IAtomContainer> moles_list = phaseItransformer(original_mole_no_h,key);
				System.out.println(key);
//				System.exit(0); 
			}
	    }
		
		
	}

}









