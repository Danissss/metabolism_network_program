package xuan.drug_porter.descriptorUtils;

/** 
 * Authors: Siyang Tian
 *
 */

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Properties;

import org.apache.commons.lang3.StringUtils;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.io.listener.PropertiesListener;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;

import xuan.drug_porter.descriptorUtils.ChemSearcher;
//import xuan.drug_porter.descriptorUtils.*;


public class SdfToSample {
	/**
	 * Create IAtomContainerSet containing all molecules in the inputFiles(either sdf or smiles).
	 * @param inputPath
	 * @return IAtomContainerSet containing all molecules with ExplicitHydrogens added.
	 * @throws Exception
	 */
	
	public IAtomContainerSet createIAtomContainerSet(String inputPath) throws Exception{
		MoleculeExplorer moleExp = new MoleculeExplorer();
		IAtomContainerSet moleculeSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
		SmilesParser sp = new SmilesParser(builder);
		String oneLine;
		//If the input is SMILEs
		if(inputPath.contains(".csv")){
			FileReader fr = new FileReader(inputPath);
			BufferedReader br = new BufferedReader(fr);
			while((oneLine = br.readLine())!=null){
				//Create IAtomContainer molecule from smiles
				IAtomContainer mol = sp.parseSmiles(oneLine);
				AtomContainerManipulator.suppressHydrogens(mol);
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
				StructureDiagramGenerator sdg = new StructureDiagramGenerator();
				sdg.setMolecule(mol);
				sdg.generateCoordinates();
				IAtomContainer layedOutMol = sdg.getMolecule();
				moleculeSet.addAtomContainer(layedOutMol);
			}

		}
		//if the input file is a sdf file
		else if(inputPath.contains(".sdf")){
			 moleculeSet = readFile(inputPath);
		}
		//If the input is a SMILES string
		else if(inputPath.contains("SMILES=")){
			String inputSMILE = inputPath.replace("SMILES=", "");
			IAtomContainer mol = sp.parseSmiles(inputSMILE);//The inputPath should be a SMILEs String in this case
			System.out.println("The input SMILES string is: " + inputSMILE);
			AtomContainerManipulator.suppressHydrogens(mol);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
			StructureDiagramGenerator sdg = new StructureDiagramGenerator();
			sdg.setMolecule(mol);
			sdg.generateCoordinates();
			IAtomContainer layedOutMol = sdg.getMolecule();
			moleculeSet.addAtomContainer(layedOutMol);
		}
		IAtomContainerSet moleculeSetWithHydrogen = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		for(int i = 0; i < moleculeSet.getAtomContainerCount(); i++){
			StructureDiagramGenerator sdg =  new StructureDiagramGenerator();
			IAtomContainer mole = moleculeSet.getAtomContainer(i);
			if(moleExp.isInvalidCandidate(mole)){
				System.out.println("Molecule with index: " + (i+1) +" is not valid, so it's treated as Non-Reactants and skipped");
				continue;
			}
			
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(mole);
			sdg.setMolecule(mole);
			sdg.generateCoordinates();
			IAtomContainer moleHyd = sdg.getMolecule();
			moleculeSetWithHydrogen.addAtomContainer(moleHyd);
			
		}
		return moleculeSetWithHydrogen;
	}
			
	/**
	 * Generate Instances for all input molecules with all features.
	 * @param set
	 * @return Instances for all input molecules with all features
	 * @throws Exception
	 */

	public Instances generateAllInstances(IAtomContainerSet set) throws Exception{
		
		 ArrayList<Attribute> atts = new ArrayList<Attribute>();
		 String moleculeFeatures = "nHBAcc\tnHBDon\tnaAromAtom\tnAtomP\tnB\tnAromBond\tnRotB\tALogP\tALogp2\tAMR\tXLogP\tMLogP\tapol\tTopoPSA\tMW\tbpol\tATSc1\tATSc2\tATSc3\tATSc4\tATSc5\tATSm1\tATSm2\tATSm3\tATSm4\tATSm5\tnAcid\tnBase\tMOMI-X\tMOMI-Y\tMOMI-Z\tMOMI-XY\tMOMI-XZ\tMOMI-YZ\tMOMI-R\tAllSurfaceArea";
		 LinkedHashMap<String, String> fpatterns = ChemSearcher.getRINFingerprintPatterns();
		 String[] labels = fpatterns.keySet().toArray(new String[fpatterns.size()]);
		 		
		 String rinFPnames = "\t" + StringUtils.join(labels,"\t");
		
		 String firstNames = moleculeFeatures+rinFPnames;
		 String[] fnames = firstNames.split("\t");
		 
		 for(int j = 0; j<fnames.length; j++){
			 fnames[j] = fnames[j].replace(",", "-");
			 Attribute Attribute = new Attribute(fnames[j]);
			 atts.add(Attribute);
		 }
		 for(int h = 0; h < 881; h++){
				
			 Attribute Attribute = new Attribute(String.format("pubchem_f%03d", h+1));
			 atts.add(Attribute);
		 }

		 for(int h = 0; h < 166; h++){
			 Attribute Attribute = new Attribute(String.format("maccs_k%03d", h+1));
			 atts.add(Attribute);
		 }
		//Add class attribute as weka request this when predicting
		Attribute classAtt = new Attribute("class"); 
		atts.add(classAtt);
		//Names have been added, now add values 
		Instances userinput = new Instances("Rel", atts, 100000);
		int length = atts.size();
		for(int idx = 0; idx<set.getAtomContainerCount();idx++){
		  System.out.println("Processing Molecue: " + (idx+1));
		  String result = generateOneinstanceFeatures(set.getAtomContainer(idx));
		  String[] temp = result.split("\t");
		  Instance sample = new DenseInstance(length); 
		  for(int vidx = 0; vidx < temp.length; vidx++){
			  Attribute att = atts.get(vidx);
			  //att.type();
			  double vle = Double.parseDouble(temp[vidx]);
			  //FastVector attributes = new FastVector();
			  sample.setValue(att, vle);
		  }
		  sample.setValue(classAtt, 0.0);
		  userinput.add(sample);
		}
		return userinput;			
	}
	
	
	/**
	 * Read the sdf file and generate a IAtomContainerSet that contains all molecules in it.
	 * @param pathToInputFile
	 * @return
	 * @throws FileNotFoundException
	 * @throws CDKException
	 */

	public IAtomContainerSet readFile(String pathToInputFile)
			throws FileNotFoundException, CDKException {
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		IteratingSDFReader sdfr = new IteratingSDFReader(new FileReader(pathToInputFile),
				bldr);
		Properties prop = new Properties();
		prop.setProperty("ForceReadAs3DCoordinates", "true");
		PropertiesListener listener = new PropertiesListener(prop);
		sdfr.addChemObjectIOListener(listener);
		sdfr.customizeJob();
		IAtomContainerSet MOLS = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		while (sdfr.hasNext())
				MOLS.addAtomContainer(sdfr.next());
		return MOLS;

	}

	/**
	 * Generate one String that contains all 2279 features for one molecule
	 * @param mole: One IAtomContainer that contains one molecule
	 * @return
	 * @throws Exception
	 */

	public String generateOneinstanceFeatures(IAtomContainer mole) throws Exception {
		StringBuffer sb = new StringBuffer();
		ChemSearcher cs = new ChemSearcher();
		PubchemFingerprinter pbf 	= new PubchemFingerprinter(SilentChemObjectBuilder.getInstance());
		MACCSFingerprinter maccs 	=  new MACCSFingerprinter(SilentChemObjectBuilder.getInstance());

		LinkedHashMap<String, String> fpatterns = cs.getRINFingerprintPatterns();
		FeatureGeneration fgen = new FeatureGeneration();
		
		IAtomContainer container = mole;
			
		IAtomContainer prepContainer = MoleculeExplorer.preprocessContainer(container);
		String[] gg = fgen.generateExtendedMolecularFeatures(prepContainer).split(",");
			
		String extendedFeatures = StringUtils.join(fgen.generateExtendedMolecularFeatures(prepContainer).split(","), "\t");
			

		ArrayList<Double> bioTransformerFingerprint_bits = cs.generateClassyfireFingerprintAsDouble(prepContainer, fpatterns).getBitValues();
		for(int x = 0; x < bioTransformerFingerprint_bits.size(); x++){
			extendedFeatures =  extendedFeatures + "\t" + String.valueOf(bioTransformerFingerprint_bits.get(x));
		}
		
			
		ArrayList<Double> fingerprint_bits = new ArrayList<Double>();
		IBitFingerprint fingerp	= pbf.getBitFingerprint(prepContainer);

		int[] onbits = fingerp.getSetbits();

		for(int kp = 0; kp < 881; kp++){
			fingerprint_bits.add(0.0);
		}
		for(int o = 0; o < onbits.length; o++){
			fingerprint_bits.set(onbits[o], 1.0);
		}
		
		extendedFeatures =  extendedFeatures + "\t" + StringUtils.join(fingerprint_bits,"\t");
			
		ArrayList<Double> maccs_fingerprint_bits = new ArrayList<Double>();
		IBitFingerprint maccs_fingerp		= maccs.getBitFingerprint(prepContainer);
			
		int[] maccs_onbits = maccs_fingerp.getSetbits();
			
		for(int kp = 0; kp < 166; kp++){
			maccs_fingerprint_bits.add(0.0);
		}
		for(int o = 0; o < maccs_onbits.length; o++){
			maccs_fingerprint_bits.set(maccs_onbits[o], 1.0);
		}			
		
		extendedFeatures =  extendedFeatures + "\t" + StringUtils.join(maccs_fingerprint_bits,"\t");
		
		String finalFeatureValues = extendedFeatures;
		String[] temp = extendedFeatures.split("\t");
		return finalFeatureValues;
	}
	
}


