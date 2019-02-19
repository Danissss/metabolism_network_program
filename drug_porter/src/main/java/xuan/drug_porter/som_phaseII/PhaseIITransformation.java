package xuan.drug_porter.som_phaseII;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.openscience.cdk.Atom;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.modeling.builder3d.ModelBuilder3D;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.query.SMARTSException;
import ambit2.smarts.query.SmartsPatternCDK;
import weka.core.Instances;
import xuan.drug_porter.som_phaseI.SomPrediction;
import xuan.drug_porter.som_phaseII.FunctionalGroup;

public class PhaseIITransformation {

static IAtom oxygen = new Atom(8);
	
	public static FunctionalGroup FG = new FunctionalGroup();
	public static SmilesGenerator smigen = new SmilesGenerator(SmiFlavor.Isomeric);
	public static String OH_SMARTS = "[OH]";
	public static String NH2_SMARTS = "N";
	
	/**
	 * 
	 * @param smiles
	 * @return
	 * @throws CDKException
	 * @throws CloneNotSupportedException
	 * @throws IOException
	 */
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
	 * @param mole
	 * @param som
	 * @param enzyme_type
	 * @return
	 * @throws Exception 
	 */
	public static ArrayList<IAtomContainer> phaseIItransformer(IAtomContainer mole, String enzyme_type) throws Exception {
		
		int num_of_atom = mole.getAtomCount();
		int num_of_bond = mole.getBondCount();
		ArrayList<IAtomContainer> biotransformered = new ArrayList<IAtomContainer>();
//		IAtom site_of_atom = mole.getAtom(som);
//		List<IAtom> connected_atom = mole.getConnectedAtomsList(site_of_atom);
//		for(int i = 0; i < connected_atom.size(); i++) {
//			
//			// rule system here: add_enzyme;
//			
//		}
//		
		
		if(enzyme_type == "UGT") {
			SmartsPatternCDK smarts = new SmartsPatternCDK();
			smarts.setSmarts(OH_SMARTS);
			// smarts.hasSMARTSPattern(mole) return number of encountered smarts
			if(smarts.hasSMARTSPattern(mole) != 0) {
				List<List<Integer>> matach_atom_ind = smarts.getUniqueMatchingAtoms(mole); // two matching group will be list.size() = 2
				for(int i = 0; i< matach_atom_ind.size(); i++) {
					List<Integer> tmp = matach_atom_ind.get(i);			// usually the first index is the som;
					int som = tmp.get(0);
					IAtomContainer new_mole = FunctionalGroup.AddPhase2Group(mole, enzyme_type, som);
					biotransformered.add(new_mole);

				}
			}
			
		}else if(enzyme_type == "SULT") {
			
			// COMT: add new sulfuate group to the OH; issue remain at which OH.
			SmartsPatternCDK smarts = new SmartsPatternCDK();
			smarts.setSmarts(OH_SMARTS);
			// smarts.hasSMARTSPattern(mole) return number of encountered smarts
			if(smarts.hasSMARTSPattern(mole) != 0) {
				List<List<Integer>> matach_atom_ind = smarts.getUniqueMatchingAtoms(mole); // two matching group will be list.size() = 2
				for(int i = 0; i< matach_atom_ind.size(); i++) {
					List<Integer> tmp = matach_atom_ind.get(i);			// usually the first index is the som;
					int som = tmp.get(0);
					IAtomContainer new_mole = FunctionalGroup.AddPhase2Group(mole, enzyme_type, som);
					biotransformered.add(new_mole);

				}
			}
			
		}
		else if(enzyme_type == "NAT") {
			
			SmartsPatternCDK smarts = new SmartsPatternCDK();
			smarts.setSmarts(NH2_SMARTS);
			// smarts.hasSMARTSPattern(mole) return number of encountered smarts
//			System.out.println(smarts.hasSMARTSPattern(mole));
			if(smarts.hasSMARTSPattern(mole) != 0) {
				List<List<Integer>> matach_atom_ind = smarts.getUniqueMatchingAtoms(mole); // two matching group will be list.size() = 2
				for(int i = 0; i< matach_atom_ind.size(); i++) {
					List<Integer> tmp = matach_atom_ind.get(i);			// usually the first index is the som;
					int som = tmp.get(0);
					
					IAtomContainer new_mole = FunctionalGroup.AddPhase2Group(mole, enzyme_type, som);
//					String output = smigen.create(new_mole);
//					System.out.println(output);
					biotransformered.add(new_mole);

				}
			}
		}
		else if(enzyme_type == "GST") {
			
		}
		else if(enzyme_type == "COMT") {
			// COMT: add new Methyl group to the OH; issue remain at which OH.
			SmartsPatternCDK smarts = new SmartsPatternCDK();
			smarts.setSmarts(OH_SMARTS);
			// smarts.hasSMARTSPattern(mole) return number of encountered smarts
			if(smarts.hasSMARTSPattern(mole) != 0) {
				List<List<Integer>> matach_atom_ind = smarts.getUniqueMatchingAtoms(mole); // two matching group will be list.size() = 2
				for(int i = 0; i< matach_atom_ind.size(); i++) {
					List<Integer> tmp = matach_atom_ind.get(i);			// usually the first index is the som;
					int som = tmp.get(0);
					IAtomContainer new_mole = FunctionalGroup.AddPhase2Group(mole, enzyme_type, som);
					biotransformered.add(new_mole);

				}
			}
		}
		else if(enzyme_type == "TAU") {
			
		}
		else if(enzyme_type == "GLYAT") {
			
		}
		
		
		return biotransformered;
	}
	
	
	
	
	public static void main(String[] args) throws Exception {
		
		String smile = "CCCCCCOC1=CC=C(N)C=C1";   // O-demethylation
//		String smile2OH = "OC1=C(O)C=C(C=C1)C(=O)\C=C\C1=CC=CC=C1";
//		[OH]
		IAtomContainer original_mole = Get3DAtomContainerFromSmiles(smile);
		IAtomContainer original_mole_no_h = AtomContainerManipulator.removeHydrogens(original_mole);
		ArrayList<IAtomContainer> transformed = phaseIItransformer(original_mole_no_h,"NAT");

		
	}
}
