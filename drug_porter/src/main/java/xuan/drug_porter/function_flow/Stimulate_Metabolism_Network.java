package xuan.drug_porter.function_flow;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
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
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.io.listener.PropertiesListener;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import weka.core.Instances;
import xuan.drug_porter.RunClassification;
import xuan.drug_porter.reactant.phase1inhibitor;
import xuan.drug_porter.reactant.phase1react;
import xuan.drug_porter.reactant.phase2react;
import xuan.drug_porter.som_phaseI.PhaseITransformation;
import xuan.drug_porter.som_phaseI.SomPrediction;
import xuan.drug_porter.som_phaseII.FunctionalGroup;



public class Stimulate_Metabolism_Network {
	
	
	
	
	public static String current_dir = System.getProperty("user.dir");
	
	/**
	 * this is for testing the system
	 * @param input
	 */
	public static void run_test(String input) {
		
		
		
		
		
		
		
	}
	
	
	
	/**
	 * this is for run real stimulation
	 * @param input
	 */
	public static void run_stimulation(String input) {
		
		
		
		
		
		
		
	}
	
	/**
	 * test for biotransformation
	 * @param input
	 * @throws Exception 
	 */
	public static void test_for_som_transformation(String input) throws Exception {
		
		
		// determine if the AtomContainer has 3D coorinator
		
//		IAtomContainer mole = get_atomcontainer_from_smiles(input);
//		boolean has_3d_con = GeometryTools.has3DCoordinates(mole);
//		
//		if (has_3d_con == false) {
//			System.out.println("Can't forward due to lack of 3D molecule!");
//			System.exit(-1);
//		}
		
		// transportation test 
		// import into blood stream from intestine
		
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
			enter_enterocyte_result.put(key, ABCC3Result.get(key));
			if(ABCC3Result.get(key).contains("non")) {
				System.out.println("the compound can't across enterocytes to enter blood stream");
			}
			
	    }
		
		
		// enter liver cell for metabolism 
		String[] enter_liver_cell_port = new String[] {"SLC22A1","SLC10A2","SLC15A1","SLCO1A2","SLCO1B1","SLCO1B3","SLCO2A1","SLCO2B1"};
		
		
		
		
		
		// cypreact
		
		
		
		// som
		
		
		// phase I transformation
		
		
		
		// phase II react
		
		
		
		// phase II transformation
		
		
		// exportation from liver to blood 
		
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
			IAtomContainer mole = PhaseITransformation.phaseItransformer(mole, som);
			metabolites_phase1.add(mole);
			
		}
		
		
			
	}
	
	
	
	/**
	 * get IAtomContainer from smiles
	 * @param smiles
	 * @return
	 * @throws CDKException
	 */
	public static IAtomContainer get_atomcontainer_from_smiles(String smiles) throws CDKException {
		
		
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
	public static String get_sdf_3D(String smiles) throws IOException {
		
		
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

	public IAtomContainer get_atomcontainer_from_sdf(String pathToInputFile) throws FileNotFoundException, CDKException {
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
		
		test_for_som_transformation("/Users/xuan/Desktop/HMDB00001.sdf");
		String sdf_file_location = get_sdf_3D("CC(=O)N(C1=C(C)C(C)=NO1)S(=O)(=O)C1=CC=C(N)C=C1");
	}

}
