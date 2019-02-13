package xuan.drug_porter.function_flow;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Properties;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.io.listener.PropertiesListener;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.modeling.builder3d.ModelBuilder3D;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import com.opencsv.CSVReader;

import weka.core.Instances;
import xuan.drug_porter.RunClassification;
import xuan.drug_porter.reactant.phase1inhibitor;
import xuan.drug_porter.reactant.phase1react;
import xuan.drug_porter.reactant.phase2react;
import xuan.drug_porter.som_phaseI.PhaseITransformation;
import xuan.drug_porter.som_phaseI.SomPrediction;
import xuan.drug_porter.som_phaseII.FunctionalGroup;
import xuan.drug_porter.som_phaseII.SomPredictionRule;



public class Stimulate_Metabolism_Network {
	
	
	
	
	public static String current_dir = System.getProperty("user.dir");
	public static SmilesGenerator smigen = new SmilesGenerator(SmiFlavor.Isomeric);


	
	/**
	 * this is for testing the system
	 * @param input
	 */
	public static void RunTest(String input) {
		
		
		
		
		
		
		
	}
	
	
	
	/**
	 * this is for run real stimulation
	 * @param input
	 */
	public static void RunStimulation(String input) {
		
		
		
		
		
		
		
	}
	
	/**
	 * test for biotransformation
	 * @param input
	 * @throws Exception 
	 */
	public static void TestForSomPrediction(String input) throws Exception {
		
		
		// determine if the AtomContainer has 3D coorinator
		IAtomContainer input_compound = GetAtomContainerFromSmiles(input);
//		IAtomContainer mole = get_atomcontainer_from_smiles(input);
//		boolean has_3d_con = GeometryTools.has3DCoordinates(mole);
//		
//		if (has_3d_con == false) {
//			System.out.println("Can't forward due to lack of 3D molecule!");
//			System.exit(-1);
//		}
		
		// transportation test 
		// import into blood stream from intestine
		
		
		// keep record of the compound 
		MetaboliteObject MO_input_compound = new MetaboliteObject(input);
		
		// keep record of all possible metabolites
		ArrayList<String> possible_metabolite_all = new ArrayList<String>();
		possible_metabolite_all.add(input);
		
		
		Instances substrate_instance = RunClassification.generate_test_instance(input, "substrate");
		Instances inhibitor_instance = RunClassification.generate_test_instance(input, "inhibitor");
		
		
		String[] blood_stream_ports = new String[] {"SLC16A1","SLC10A2","SLC15A1","SLCO1A2","SLCO1B1","SLCO1B3","SLCO2A1","SLCO2B1"}; // enter enterocyte cell: cell line on intestine wall
		HashMap<String,String> enter_enterocyte_result = new HashMap<String,String>();
		int enter_enterocyte_number = 0;
		for(int i = 0; i < blood_stream_ports.length; i++) {
			HashMap<String,String> single_protein_result = RunClassification.run_classifier(substrate_instance, blood_stream_ports[i]);
			for (String key : single_protein_result.keySet()) {
				enter_enterocyte_result.put(key, single_protein_result.get(key));
				if(!single_protein_result.get(key).contains("non")) {
					enter_enterocyte_number++;
				}
				
		    }
		}
		
		if(enter_enterocyte_number > 1) {
			System.out.println("the compound enter enterocytes");
		}else {
			System.out.println("the compound can't enter enterocytes");
		}
		
		// enter blood stream 
		HashMap<String,String> ABCC3Result = RunClassification.run_classifier(substrate_instance, "ABCC3"); // only ABCC3 facilicate export from enterocyte to blood stream
		for (String key : ABCC3Result.keySet()) {
			if(ABCC3Result.get(key).contains("non")) {
				System.out.println("the compound can't across enterocytes to enter blood stream");
			}
			
	    }
		
		
		// enter liver cell for metabolism 
		
		String[] enter_liver_cell_port = new String[] {"SLC22A1","SLC10A2","SLC15A1","SLCO1A2","SLCO1B1","SLCO1B3","SLCO2A1","SLCO2B1"};
		HashMap<String,String> enter_liver_result = new HashMap<String,String>();
		int enter_liver_number = 0;
		for (int i = 0; i < enter_liver_cell_port.length; i++) {
			HashMap<String,String> tmp_enter_liver_result = RunClassification.run_classifier(substrate_instance, enter_liver_cell_port[i]); // only ABCC3 facilicate export from enterocyte to blood stream
			for (String key : tmp_enter_liver_result.keySet()) {
				enter_liver_result.put(key, tmp_enter_liver_result.get(key));
				if(!tmp_enter_liver_result.get(key).contains("non")) {
					enter_liver_number++;
				}
				
		    }
		}
		
		if(enter_liver_number > 1) {
			System.out.println("the compound enter liver cell for potential metabolism");
		}else {
			System.out.println("the compound can't enter liver cell");
			// if can't enter the liver cell, then create the branch that determine if the compound will be execrted safely
		}
		
		
		// if the compound can enter cypreact
		String[] cyp_enzyme_list = new String[] {"CYP1A2","CYP2B6","CYP2A6","CYP2C8","CYP2C9","CYP2C19","CYP2D6","CYP2E1","CYP3A4"};
		ArrayList<HashMap<String,String>> cypreact_result = phase1react.makePrediction("", String.format("SMILES=%?",input), "CYP1A2,CYP2B6,CYP2A6,CYP2C8,CYP2C9,CYP2C19,CYP2D6,CYP2E1,CYP3A4");
		for(int i = 0; i< cypreact_result.size(); i++) {
			HashMap<String,String> single_result = cypreact_result.get(i);
			for (String key : single_result.keySet()) {
//				enter_enterocyte_result.put(key, tmp_enter_liver_result.get(key));
//				if(!tmp_enter_liver_result.get(key).contains("non")) {
//					enter_liver_number++;
//				}
				System.out.println(key+" : " + single_result.get(key));
				
				
				
				
//				
		    }
		}
		
		
		// som and phase I transformation
		Instances som_instance = SomPrediction.create_test_instance(input);
		
		// for substrate of CYP, run the som prediction
		String[] substrate_canadidates = new String[9];
		ArrayList<IAtomContainer> result_of_phaseI_substrate = new ArrayList<IAtomContainer>();
		ArrayList<String> result_of_phaseI_substrate_smiles = new ArrayList<String>();
		for(int i = 0; i < substrate_canadidates.length; i++) {
			HashMap<Integer,String> som_phase1_result = SomPrediction.runSomClassifier(som_instance,substrate_canadidates[i]);
			// this is map of som for each cyp
			// create the new mol at this place and append to result_of_phaseI_substrate
			for (Integer key : som_phase1_result.keySet()) {
				if(som_phase1_result.get(key) == "Yes") {
					IAtomContainer canadidate_substrate = PhaseITransformation.phaseItransformer(input_compound, key);
					result_of_phaseI_substrate.add(canadidate_substrate);
					
					
					// append the metabolites of smiles to keep record
					String tmp_smiles = smigen.create(canadidate_substrate);
					result_of_phaseI_substrate_smiles.add(tmp_smiles);
					possible_metabolite_all.add(tmp_smiles);
				}
				
		    }
		}
		
		
		
		//String key_string = String.format("%s_substrate", enzyme_name);
		
		//if(result == 0.0) {
		//	classified_result.put(key_string, "non-substrate");
		//}else {
		//	classified_result.put(key_string, "substrate");
		//}
		
		// phase II react
		// run Original compound as well as compound from phaseI biotransformation
		ArrayList<HashMap<String,String>> phaseII_result = new ArrayList<HashMap<String,String>>();
		HashMap<String,String> original_compound_phaseII_result = phase2react.RunClassifierForAll(substrate_instance);
		
//		phaseII_result.add(phase2react.RunClassifierForAll(substrate_instance);
		
		for(int i = 0; i< result_of_phaseI_substrate_smiles.size(); i++) {
			Instances temp_substrate_instance = RunClassification.generate_test_instance(result_of_phaseI_substrate_smiles.get(i), "substrate");
			HashMap<String,String> temp_substrate_result_phaseII = phase2react.RunClassifierForAll(temp_substrate_instance);
			phaseII_result.add(temp_substrate_result_phaseII);
			
		}
		
		
		// phase II transformation: for different enzyme, do different function group adding;
		// prediction for input compound 
		for(String key: original_compound_phaseII_result.keySet()) {
			if(original_compound_phaseII_result.get(key) == "substrate") {
				// do transformation 
				// IAtomContainer AddPhase2Group(IAtomContainer mole,  String enzyme_name,int som_site)
				String temp_enzyme = key.split("_")[0];
				ArrayList<Integer> som_phaseII = SomPredictionRule.GetSomPhaseII(input_compound);
				for(int i = 0; i<som_phaseII.size(); i++) {
					IAtomContainer tmp_container = FunctionalGroup.AddPhase2Group(input_compound,temp_enzyme,som_phaseII.get(i));
					String tmp_smiles = smigen.create(tmp_container);
					possible_metabolite_all.add(tmp_smiles);
				}
				
			}
		}
		
		// prediction for phaseI metabolites
		for(int i = 0; i < phaseII_result.size(); i++) {
			
			HashMap<String,String> phaseI_metabolite_result = phaseII_result.get(i);
			for(String key: phaseI_metabolite_result.keySet()) {
				if(phaseI_metabolite_result.get(key) == "substrate") {
					// do transformation 
					// IAtomContainer AddPhase2Group(IAtomContainer mole,  String enzyme_name,int som_site)
					String temp_enzyme = key.split("_")[0];
					ArrayList<Integer> som_phaseII = SomPredictionRule.GetSomPhaseII(input_compound);
					for(int m = 0; m<som_phaseII.size(); m++) {
						IAtomContainer tmp_container = FunctionalGroup.AddPhase2Group(input_compound,temp_enzyme,som_phaseII.get(m));
						
						String tmp_smiles = smigen.create(tmp_container);
						possible_metabolite_all.add(tmp_smiles);
					}
				}
			}
		}
		
		
		
		
		// exportation from liver to blood == bioactive compound
		String[] export_from_liver_to_blood = new String[] {};
		for(int i = 0; i < possible_metabolite_all.size(); i++) {
			Instances all_metabolites_instance = RunClassification.generate_test_instance(possible_metabolite_all.get(i), "substrate");
			
		}
		
		
		
		//blood to urine 
		
		ArrayList<IAtomContainer> metabolites_phase1 = new ArrayList<IAtomContainer>();
		
		String smiles = input;
//		String smiles = String.format("SMILES=%s", input);
		phase1react phase1sub = new phase1react();
		ArrayList<HashMap<String,String>> result = phase1sub.makePrediction("", input, "CYP1A2");
		for(int i = 0; i < result.size(); i++) {
			String result1 = result.get(i).get("1A2");
		}
		
		SomPrediction spd = new SomPrediction();
		Instances instance = spd.create_test_instance(smiles);
		HashMap<Integer,String> som_classifier = spd.runSomClassifier(instance, "CYP1A2");
		ArrayList<Integer> som_phase1 = new ArrayList<Integer>();
		for(Integer key: som_classifier.keySet()) {
			som_phase1.add(key);
			
		}
		
		
		for(int i = 0; i < som_phase1.size(); i++) {
			int som = som_phase1.get(i);
			IAtomContainer mole = PhaseITransformation.phaseItransformer(input_compound, som);
			metabolites_phase1.add(mole);
			
		}
		
		
			
	}
	
	
	
	/**
	 * get IAtomContainer from smiles
	 * @param smiles
	 * @return
	 * @throws CDKException
	 */
	public static IAtomContainer GetAtomContainerFromSmiles(String smiles) throws CDKException {
		
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
	    IAtomContainer mole = builder.newInstance(IAtomContainer.class);
		SmilesParser temp_smiles = new SmilesParser(builder);
	 	IAtomContainer atom_container = temp_smiles.parseSmiles(smiles);
	 	AtomContainerManipulator.suppressHydrogens(atom_container);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(atom_container);
	
	 	StructureDiagramGenerator sdg = new StructureDiagramGenerator();
		sdg.setMolecule(atom_container);
		sdg.generateCoordinates();
		mole = sdg.getMolecule();
		
		return mole;
	}

	/**
	 * get 3d sdf based on the web server
	 * @param smiles
	 * @return
	 * @throws IOException
	 */
	public static String GetSdf3D(String smiles) throws IOException {
		
		
		String smiles_name = smiles.replace("/", "_");
		String file_location = String.format("%s/tmp_file/%s.sdf",current_dir,smiles_name);
		String stripped_smiles = smiles.replace("#", "%23");
		
		
		String url = String.format("https://cactus.nci.nih.gov/chemical/structure/%s/sdf", stripped_smiles);
		
		URL urls = new URL(url);
        URLConnection conn = urls.openConnection();
        BufferedReader in = new BufferedReader(new InputStreamReader(conn.getInputStream()));
        Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file_location), "utf-8"));
     
        String inputLine;

        while ((inputLine = in.readLine()) != null) { 
        		writer.write(inputLine);
        }
        in.close();
		writer.close();
		
		
		return file_location;
	}
	
	
	/**
	 * Read the sdf file and generate a IAtomContainerSet that contains all molecules in it.
	 * @param pathToInputFile
	 * @return
	 * @throws FileNotFoundException
	 * @throws CDKException
	 */

	public IAtomContainer GetAtomContainerFromSdf(String pathToInputFile) throws FileNotFoundException, CDKException {
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
		IAtomContainer moles = MOLS.getAtomContainer(0);
		
		
		return moles;

	}
	
	/**
	 * generate 3D atomcontainer from 2d sdf
	 * @param pathToInputFile
	 * @return
	 * @throws CDKException
	 * @throws IOException 
	 * @throws CloneNotSupportedException 
	 */
	public IAtomContainer Get3DAtomContainerFromSdf(String pathToInputFile) throws CDKException, CloneNotSupportedException, IOException {
		IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
		IteratingSDFReader sdfr = new IteratingSDFReader(new FileReader(pathToInputFile),
				builder);
		Properties prop = new Properties();
		prop.setProperty("ForceReadAs3DCoordinates", "true");
		PropertiesListener listener = new PropertiesListener(prop);
		sdfr.addChemObjectIOListener(listener);
		sdfr.customizeJob();
		IAtomContainerSet MOLS = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		while (sdfr.hasNext())
				MOLS.addAtomContainer(sdfr.next());
		IAtomContainer moles = MOLS.getAtomContainer(0);
		
		ModelBuilder3D mb3d = ModelBuilder3D.getInstance(builder);
		IAtomContainer molecule_3D = mb3d.generate3DCoordinates(moles, false);
		
		return molecule_3D;

	}
	
	
	/**
	 * helper function that parse the input: e.g. convert smiles to 3D container; check if the sdf is 3d and then convert it to 3d;
	 * return the IAtomContainer object
	 * @param input
	 * @return
	 */
	public IAtomContainer ParseInputToAtomContainer(String input) {
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
	    IAtomContainer mole = builder.newInstance(IAtomContainer.class);
	    
	    
	    
	    
	    
	    
	    return mole;
		
		
	}
	
	/**
	 * get 3d IAtomContainer from smiles
	 * @param smiles
	 * @return
	 * @throws CDKException
	 * @throws IOException 
	 * @throws CloneNotSupportedException 
	 */
	public static IAtomContainer Get3DAtomContainerFromSmiles(String smiles) throws CDKException, CloneNotSupportedException, IOException {
		
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
	    IAtomContainer mole = builder.newInstance(IAtomContainer.class);
		SmilesParser temp_smiles = new SmilesParser(builder);
	 	IAtomContainer atom_container = temp_smiles.parseSmiles(smiles);
	 	AtomContainerManipulator.suppressHydrogens(atom_container);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(atom_container);
	
	 	StructureDiagramGenerator sdg = new StructureDiagramGenerator();
		sdg.setMolecule(atom_container);
		sdg.generateCoordinates();
		mole = sdg.getMolecule();

		ModelBuilder3D mb3d = ModelBuilder3D.getInstance(builder);
		IAtomContainer molecule_3D = mb3d.generate3DCoordinates(mole, false);
		
		return molecule_3D;
	}
	
	
	
	/**
	 * main function
	 * @param args
	 * @throws Exception
	 */
	public static void main(String args[]) throws Exception {
		
		// convert smiles to 3d sdf
		// transporter
		// if pass; keep going; else: keep going
		// phaseIinhibitor -> generate information
		// phaseIreact -> phaseIsom -> transformation -> new molecule (phaseImetabolites) -> test for phaseIinhibitor -> 
		// phaseIIreact -> phaseIItransformation
		// phaseImetabolites -> phaseIIreact -> phaseIItransformation
		// keep all new compounds (metabolites) -> phaseI inhibitor + phaseII inhibitor + transporter exertion and inhibitation
		// report
		
		TestForSomPrediction("/Users/xuan/Desktop/HMDB00001.sdf");
//		String sdf_file_location = GetSdf3D("CC(=O)N(C1=C(C)C(C)=NO1)S(=O)(=O)C1=CC=C(N)C=C1");
		
//		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
//		ModelBuilder3D mb3d = ModelBuilder3D.getInstance(builder);
//		IAtomContainer molecule_2D = GetAtomContainerFromSmiles("CC(=O)N(C1=C(C)C(C)=NO1)S(=O)(=O)C1=CC=C(N)C=C1");
//		IAtomContainer molecule_3D = mb3d.generate3DCoordinates(molecule_2D, false);
//		
//		String tempFile = current_dir+"/tmp_file/3d_molecule.sdf";
//	 	SDFWriter sdw  = new SDFWriter(new FileWriter(tempFile));
//	 	sdw.write(molecule_3D);
//	 	sdw.close();
	}

}
