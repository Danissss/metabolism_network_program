package xuan.drug_porter.som_phaseII;

import java.util.HashMap;
import java.util.List;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.Element;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.interfaces.*;

import xuan.drug_porter.descriptorUtils.SdfToSample;

public class FunctionalGroup {
	
	
	static IAtom carbon = new Atom(6);
	static IAtom nitrogen = new Atom(7);
	static IAtom oxygen = new Atom(8);
	static IAtom gold  = new Atom(79); // this is heavy metal for indexing
	
	private static String SULT = "OS(=O)=O";	  // 1
	private static String UGT  = "O[C@H]1CO[C@@H]([C@@H](O)[C@@H]1O)C(O)=O";  // 2
	private static String NAT  = "CC=O";   // 1
	private static String GST  = "N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(O)=O)C(O)=O"; //9
	private static String COMT = "C"; // 0
	private static String TAU  = "NCCS(O)(=O)=O"; //0
	private static String GLYAT = "NCC(O)=O"; // 0
	
	
	static SdfToSample sts  = new SdfToSample();
	static HashMap<String, Integer> som_map_phaseII = MapSOM();
	
	/**
	 * 
	 * @param mole
	 * @param enzyme_name
	 * @param som_site
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainer AddPhase2Group(IAtomContainer mole,  String enzyme_name,int som_site) throws Exception {
	
		int start_index = new Integer(0);
		String smiles_string = new String();
		IAtomContainer enzyme_mole = new AtomContainer();
		SmilesGenerator smigen = new SmilesGenerator(SmiFlavor.Isomeric);
		
		if(enzyme_name == "SULT") {
			start_index = som_map_phaseII.get(SULT);
			smiles_string = SULT;
		}else if(enzyme_name == "UGT"){
			start_index = som_map_phaseII.get(UGT);
			smiles_string = UGT;
		}
		else if(enzyme_name == "NAT"){
			start_index = som_map_phaseII.get(NAT);
			smiles_string = NAT;
		}
		else if(enzyme_name == "GST"){
			start_index = som_map_phaseII.get(GST);
			smiles_string = GST;
		}
		else if(enzyme_name == "COMT"){
			start_index = som_map_phaseII.get(COMT);
			smiles_string = COMT;
		}
		else if(enzyme_name == "TAU"){
			start_index = som_map_phaseII.get(TAU);
			smiles_string = TAU;
		}
		else if(enzyme_name == "GLYAT"){
			start_index = som_map_phaseII.get(GLYAT);
			smiles_string = GLYAT;
		}
		else {
			System.out.println("Unkown enzyme type!");
			System.exit(-1);
		}
		
		IAtomContainerSet mole_set = sts.createIAtomContainerSet("SMILES="+smiles_string);
		enzyme_mole = mole_set.getAtomContainer(0); // with hydrogen
		IAtomContainer enzyme_mole_no_h = AtomContainerManipulator.removeHydrogens(enzyme_mole);
		
		
		
		enzyme_mole_no_h.getAtom(start_index).setCharge(1.0);
		mole.getAtom(som_site).setCharge(2.0);
		
		
		mole.add(enzyme_mole_no_h);
		
		int index_1 = new Integer(0);
		int index_2 = new Integer(0);
		for (int k = 0; k < mole.getAtomCount(); k++){
			IAtom atoms = mole.getAtom(k);
			if(atoms.getCharge() != null) {
				
			
			if(atoms.getCharge() == 1.0) {
				index_1 = k;
				atoms.setCharge(null);
			}
			else if (atoms.getCharge() == 2.0) {
				index_2 = k;
				atoms.setCharge(null);
			}
			}
		}
		
		// add bond / make connection
		mole.addBond(index_1,index_2,IBond.Order.SINGLE);
		
		
		// remove the charges 
		for (int k = 0; k < mole.getAtomCount(); k++){
			IAtom atoms = mole.getAtom(k);
			if(atoms.getCharge() != null) {
//				System.out.println(atoms.getCharge());
				atoms.setCharge(null);
//				System.out.println(atoms.getCharge());
			}
		}
		
		
		
		return mole;
	}
	
	
	
	private static HashMap<String, Integer> MapSOM() {
		HashMap<String, Integer> myMap = new HashMap<String, Integer>();
	    myMap.put("SULT", 1);
	    myMap.put("UGT", 2);
	    myMap.put("NAT", 1);
	    myMap.put("GST", 9);
	    myMap.put("COMT", 0);
	    myMap.put("TAU", 0);
	    myMap.put("GLYAT", 0);
	    
	    
	    return myMap;
	}
	
	
	
	public static void main(String[] args) throws Exception {
		
		IAtomContainerSet mole_set = sts.createIAtomContainerSet("SMILES=CC1=CC=CC(NC2=CC=CC=C2C(O)=O)=C1C");
		IAtomContainer mole = mole_set.getAtomContainer(0); // with hydrogen
		IAtomContainer mole_no_h = AtomContainerManipulator.removeHydrogens(mole);

		IAtomContainer ugt = AddPhase2Group(mole_no_h,"UGT",14); // 16 but minus -1  = 15
		
		SmilesGenerator smigen = new SmilesGenerator(SmiFlavor.Isomeric);
		String smi = smigen.create(ugt);

		
	}
}
