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
	 * refactor: taking input to determine whether it can enter liver or not
	 * @param input
	 * @return
	 * @throws Exception
	 */
	public static MetaboliteObject EnterLiver(String input) throws Exception{
		
		MetaboliteObject MO_input_compound = new MetaboliteObject(input,input);
		
		
		// keep record of all possible metabolites
		ArrayList<String> possible_metabolite_all = new ArrayList<String>();
		possible_metabolite_all.add(input);
		
		
		Instances substrate_instance = RunClassification.generate_test_instance(input, "substrate");
		
		
		String[] blood_stream_ports = new String[] {"SLC16A1","SLC10A2","SLC15A1","SLCO1A2","SLCO1B1","SLCO1B3","SLCO2A1","SLCO2B1"}; // enter enterocyte cell: cell line on intestine wall
		HashMap<String,String> enter_enterocyte_result = new HashMap<String,String>();

		for(int i = 0; i < blood_stream_ports.length; i++) {
			HashMap<String,String> single_protein_result = RunClassification.run_classifier(substrate_instance, blood_stream_ports[i],"substrate");
			for (String key : single_protein_result.keySet()) {
				enter_enterocyte_result.put(key, single_protein_result.get(key));
				
				if(!single_protein_result.get(key).contains("non")) {
					MO_input_compound.AddTargetAsSubstrate(blood_stream_ports[i]);
				}
				
		    }
		}
		
		if(MO_input_compound.GetAllTargetSubstrate().size() > 0) {
			System.out.println("Input compound enter enterocytes");
		}
		else {
			System.out.println("Input compound can't enter enterocytes");
		}
		
		
		
		// enter blood stream 
		HashMap<String,String> ABCC3Result = RunClassification.run_classifier(substrate_instance, "ABCC3","substrate"); // only ABCC3 facilicate export from enterocyte to blood stream
		for (String key : ABCC3Result.keySet()) {
			if(ABCC3Result.get(key).contains("non")) {
				System.out.println("Input compound can't across enterocytes to enter blood stream");
			}else {
				MO_input_compound.AddTargetAsSubstrate("ABCC3");
				System.out.println("Input compound across enterocytes to enter blood stream");
			}
			
	    }
		
		
		return MO_input_compound;
		
		
	}
	
	
	
	/**
	 * determine which transporter will be inhibited by the input molecule 
	 * if the input molecule inhibit all transporter in liver; then no metabolite may go out
	 * 
	 * @param input
	 * @return
	 * @throws Exception
	 */
	public static void InhibitLiverTransporter(MetaboliteObject input) throws Exception{
		
		
		Instances inhibitor_instance = RunClassification.generate_test_instance(input.SMILES, "inhibitor");
		String[] liver_transporter = new String[] {"SLC22A1","SLC10A2","SLC15A1","SLCO1A2","SLCO1B1","SLCO1B3","SLCO2A1","SLCO2B1","ABCC3","ABCC4","ABCC6","ABCB1","ABCG2","ABCB11","ABCC2","SLC47A1"}; 
		for (int i = 0; i<liver_transporter.length; i++) {
			HashMap<String,String> single_result = RunClassification.run_classifier(inhibitor_instance, liver_transporter[i],"inhibitor");
			for (String key : single_result.keySet()) {
				single_result.put(key, single_result.get(key));
				
				if(!single_result.get(key).contains("non")) {
					input.AddTargetAsInhibitor((liver_transporter[i]));
				}
				
		    }
			
		}
		
	}
	
	
	
	
	/**
	 * determine which enzyme will be inhibited by input molecule
	 * if all enzyme is inhibited; then no metabolism; 
	 * caution: enzyme may cause conflict 
	 * TODO: add phase II inhibition
	 * @param input
	 * @return
	 * @throws Exception
	 */
	public static void InhibitLiverEnzyme(MetaboliteObject input) throws Exception {
		
		Instances single_inhibitor_instance = RunClassification.generate_test_instance(input.SMILES, "inhibitor");
		String[] cyp_enzyme_list = new String[] {"CYP1A2","CYP2B6","CYP2A6","CYP2C8","CYP2C9","CYP2C19","CYP2D6","CYP2E1","CYP3A4"};
		
		for (int i = 0; i < cyp_enzyme_list.length; i++) {
			HashMap<String,String> single_result = phase1inhibitor.run_classifier(single_inhibitor_instance, cyp_enzyme_list[i]);
			for (String key : single_result.keySet()) {
				single_result.put(key, single_result.get(key));
				if(!single_result.get(key).contains("non")) {
					input.AddTargetAsInhibitor((cyp_enzyme_list[i]));
				}
				
		    }
		}
		
		
	}
	
	
	
	
	/**
	 * 
	 * @param input
	 * @throws Exception
	 */
	public static ArrayList<MetaboliteObject> LiverMetabolismPhaseI(MetaboliteObject input) throws Exception {
		
//		ArrayList<MetaboliteObject> phaseImetabolites = new ArrayList<MetaboliteObject>();
		
		
		String[] cyp_enzyme_list = new String[] {"CYP1A2","CYP2B6","CYP2A6","CYP2C8","CYP2C9","CYP2C19","CYP2D6","CYP2E1","CYP3A4"};
		ArrayList<HashMap<String,String>> cypreact_result = phase1react.makePrediction("", String.format("SMILES=%?",input), "CYP1A2,CYP2B6,CYP2A6,CYP2C8,CYP2C9,CYP2C19,CYP2D6,CYP2E1,CYP3A4");
		for(int i = 0; i< cypreact_result.size(); i++) {
			HashMap<String,String> single_result = cypreact_result.get(i);
			for (String key : single_result.keySet()) {
				if(single_result.get(key) == "R") {
					// TODO: need to change the return index
					input.AddTargetAsSubstrate(cyp_enzyme_list[i]);
				}
		    }
		}
		
		
		// som and phase I transformation
		// for substrate of CYP, run the som prediction
		// find which one is reactant, do the som prediction, do the phaseI transformation
		ArrayList<MetaboliteObject> phaseImetabolites = new ArrayList<MetaboliteObject>();
		Instances som_instance = SomPrediction.create_test_instance(input.SMILES);
		for(int i = 0; i < cyp_enzyme_list.length; i++) {
			if(input.GetAllTargetSubstrate().contains(cyp_enzyme_list[i])) {
				// it is substrate for this cyp
				HashMap<Integer,String> som_phase1_result = SomPrediction.runSomClassifier(som_instance,cyp_enzyme_list[i]);
				for (Integer key : som_phase1_result.keySet()) {
					if(som_phase1_result.get(key) == "Yes") {
						IAtomContainer mol = Get3DAtomContainerFromSmiles(input.SMILES);
						ArrayList<IAtomContainer> canadidate_substrate = PhaseITransformation.phaseItransformer(mol, key);
						for(int k = 0; k < canadidate_substrate.size(); k++) {
							String tmp_smiles = smigen.create(canadidate_substrate.get(k));
							MetaboliteObject tmp_obj = new MetaboliteObject(tmp_smiles,input.SMILES);
										
							// change the objects for parents
							phaseImetabolites.add(tmp_obj);
							input.AddChild(tmp_smiles);
							
						}
					}
					
			    }
			}
		}
		
		return phaseImetabolites;
		
	}
	
	
	
	/**
	 * phase II react
	 * and phase II biotransformation
	 * run Original compound as well as compound from phaseI biotransformation
	 * @param input
	 * @return
	 * @throws Exception 
	 */
	public static ArrayList<MetaboliteObject> LiverMetabolismPhaseII(ArrayList<MetaboliteObject> input) throws Exception{
		// phase II react
		// and phase II biotransformation
		// run Original compound as well as compound from phaseI biotransformation
		ArrayList<MetaboliteObject> phaseIImetabolites = new ArrayList<MetaboliteObject>();
		
		// phaseII_result.add(phase2react.RunClassifierForAll(substrate_instance);
		// each MetaboliteObject smiles will undergo phaseII_model_name
		String[] phaseII_model_name = new String[] {"UGT","SULT","NAT","GST","COMT"};
		for(int i = 0; i< input.size(); i++) {
			String current_metabolite_smiles = input.get(i).SMILES;
			Instances temp_substrate_instance = RunClassification.generate_test_instance(current_metabolite_smiles, "substrate");
			IAtomContainer phaseImetabolite = GetAtomContainerFromSmiles(input.get(i).SMILES);
			
			// get the result as <UGT, non-substrate/substrate>
			for(int k = 0; k < phaseII_model_name.length; k++) {
				HashMap<String,String> temp_substrate_result_phaseII = phase2react.RunClassifier(temp_substrate_instance,phaseII_model_name[k]);
				for(String key: temp_substrate_result_phaseII.keySet()) {
					if(temp_substrate_result_phaseII.get(key) == "substrate") {
						// it is the substrate of the current enzyme 
					    System.out.println(String.format("% is substrate of %s.", input.get(i).SMILES,key));
						input.get(i).AddTargetAsSubstrate(phaseII_model_name[k]);
						IAtomContainer mole = GetAtomContainerFromSmiles(input.get(i).SMILES); 			// 2D for now. For future if SOM predictor here change it to 3D
						
						// do transformation
						ArrayList<Integer> som_phaseII = SomPredictionRule.GetSomPhaseII(mole);
						for(int m = 0; m < som_phaseII.size(); m++) {
							IAtomContainer tmp_container = FunctionalGroup.AddPhase2Group(phaseImetabolite,phaseII_model_name[k],som_phaseII.get(m));
							String tmp_smiles = smigen.create(tmp_container);
							MetaboliteObject tmp_phaseII = new MetaboliteObject(tmp_smiles,current_metabolite_smiles);
							phaseIImetabolites.add(tmp_phaseII);
							
							// add the child for current phaseImetabolite object
							input.get(i).AddChild(tmp_smiles);
						}				
					}
				}
			}
		}
		
		return phaseIImetabolites;
		
		
		
	}
	
	/**
	 * execration from liver to blood == bioactive compound
	 * only care about the bioactive stuff; metabolites go to bile is bioinactive
	 * determine which metabolite will be efflux and become bioactive compound
	 * tranpsorter to export to blood (major) = new String[] {"ABCC3","ABCC4","ABCC6"};	
	 * @param input
	 * @return
	 * @throws Exception 
	 */
	public static ArrayList<MetaboliteObject> BioactiveMetabolites(ArrayList<MetaboliteObject> input) throws Exception{
		
		
		// or execration from liver to bile == bio-inactive compound
		String[] export_from_liver_to_blood = new String[] {"ABCC3","ABCC4","ABCC6"};
		ArrayList<MetaboliteObject> export_from_liver_to_blood_result = new ArrayList<MetaboliteObject>();  				
		for(int i = 0; i < input.size(); i++) {
			String current_metabolite_smiles = input.get(i).GetSMILES();
			Instances single_metabolites_instance = RunClassification.generate_test_instance(current_metabolite_smiles, "substrate");
			// for each metabolites, test for each transporters
			int bioactive_flag = 0;
			
			for (int k = 0; k < export_from_liver_to_blood.length; k++) {
				HashMap<String,String> single_export_result = RunClassification.run_classifier(single_metabolites_instance, export_from_liver_to_blood[i],"substrate");
				for (String key : single_export_result.keySet()) {
					if(single_export_result.get(key) == "substrate") {
						System.out.println(String.format("%s export from liver to blood by %s (bioactive).", current_metabolite_smiles, export_from_liver_to_blood[i]));
						input.get(i).BioActive(true);
						
						// add the transporter as substrate
 						if(!input.get(i).GetAllTargetSubstrate().contains(export_from_liver_to_blood[i])) {
 							input.get(i).AddChild(export_from_liver_to_blood[i]);
						}
 						
 						bioactive_flag = 1;
						
					}
					else {
						System.out.println(String.format("%s doesn't export from liver to blood (bioactive).", current_metabolite_smiles));
					}
			    }
			}
			
			
			if (bioactive_flag == 1) {
				export_from_liver_to_blood_result.add(input.get(i));
				
			}
		}
		
		
		return export_from_liver_to_blood_result;
		
		
	}
	
	
	
	/**
	 * exportation from liver to bile == non-bioactive compound
	 * determine which metabolite will be efflux from liver  
	 * transporter for exporting to bile = {"ABCB1","ABCG2","ABCB11","ABCC2","SLC47A1"};    
	 * @param input
	 * @return
	 * @throws Exception 
	 */
	public static ArrayList<MetaboliteObject> EfflexMetabolite(ArrayList<MetaboliteObject> input) throws Exception{

		
		String[] export_from_liver_to_bile = new String[] {"ABCB1","ABCG2","ABCB11","ABCC2","SLC47A1"};    				
		ArrayList<MetaboliteObject> export_from_liver_to_bile_result = new ArrayList<MetaboliteObject>();  									
		for(int i = 0; i < input.size(); i++) {
			int export_flag = 0;

			
			String current_metabolite_smiles = input.get(i).GetSMILES();
			Instances single_metabolites_instance = RunClassification.generate_test_instance(current_metabolite_smiles, "substrate");
			// for each metabolites, test for each transporters
			for (int k = 0; k < export_from_liver_to_bile.length; k++) {
				HashMap<String,String> single_export_result = RunClassification.run_classifier(single_metabolites_instance, export_from_liver_to_bile[i],"sub");
				for (String key : single_export_result.keySet()) {
					if(single_export_result.get(key) == "substrate") {
						System.out.println(String.format("%s export from liver to bile by %s (bio-inactive).", current_metabolite_smiles,export_from_liver_to_bile[i]));
						
						// add the transporter as substrate
						if(!input.get(i).GetAllTargetSubstrate().contains(export_from_liver_to_bile[i])) {
							input.get(i).AddChild(export_from_liver_to_bile[i]);
						}
						export_flag = 1;
					
					}
					
					else {
						System.out.println(String.format("%s doesn't export from liver to blood (bioactive).", current_metabolite_smiles));
					}
			    }
			}
			
			
			
			
			
			if (export_flag == 1) {
				export_from_liver_to_bile_result.add(input.get(i));
			}
			
		}
		
		
		//String[] exports_ = new String[] {"ABCC3","ABCC4","ABCC6","ABCB1","ABCG2","ABCB11","ABCC2","SLC47A1"};	
		// for test if all metabolites can be exported, test if they cntain one of these {"ABCC3","ABCC4","ABCC6","ABCB1","ABCG2","ABCB11","ABCC2","SLC47A1"} transporters.
		return export_from_liver_to_bile_result;
		
	}
	
	
	
	
	/**
	 * determine which metabolite will be efflux from urine 
	 * only accept the list with bioactive compounds
	 * return the metabolite that will be exported from urine; some other may re-enter liver for bile exertion
	 * @param input
	 * @throws Exception 
	 */
	public static ArrayList<MetaboliteObject> UrineMetabolite(ArrayList<MetaboliteObject> input) throws Exception{
		
		
		
		// Blood to urine 
		// metabolites will go from blood to kidney and then urine
		// for those couldn't enter kidney, may also cause trouble
		String[] enter_kidney_tubules = new String[] {"SLCO4C1","SLC22A6","SLC22A7","SLC22A8"};
		ArrayList<MetaboliteObject> enter_kidney_tubules_result = new ArrayList<MetaboliteObject>();
		for(int i = 0; i <input.size(); i++) {
			Instances all_metabolites_instance = RunClassification.generate_test_instance(input.get(i).SMILES, "substrate");
			int enter_kidney_flag = 0;
			for(int k = 0; k < enter_kidney_tubules.length; k++) {
				HashMap<String,String> single_import_result = RunClassification.run_classifier(all_metabolites_instance, enter_kidney_tubules[i],"substrate");
				for (String key : single_import_result.keySet()) {
					if(single_import_result.get(key) == "substrate") {
						System.out.println(String.format("%s import to kidney tubules by %s.", input.get(i).SMILES,enter_kidney_tubules[i]));
						enter_kidney_tubules_result.add(input.get(i));
						enter_kidney_flag = 1;
					}
					else {
						System.out.println(String.format("%s doesn't export from liver to blood (bioactive).", input.get(i)));
					}
			    }
			}
		}
		
		
		// kidney tubules to urine
		// don't care about reabsorption
		ArrayList<MetaboliteObject> UrineMetabolite = new ArrayList<MetaboliteObject>();
		if (enter_kidney_tubules_result.size() == 0) {
			System.out.println("Nothing enter kidney.");

		}else {
			String[] enter_urine = new String[] {"SLC22A11","SLC22A12","ABCC2","ABCC4","SLC47A1","SLC47A2","ABCB1","SLC22A4","SLC22A5"};
			for(int i = 0; i <enter_kidney_tubules_result.size(); i++) {
				
				
				Instances all_metabolites_instance = RunClassification.generate_test_instance(enter_kidney_tubules_result.get(i).SMILES, "substrate");
				
				int export_to_urine_flag = 0;
				for(int k = 0; k < enter_urine.length; k++) {
					HashMap<String,String> single_import_result = RunClassification.run_classifier(all_metabolites_instance, enter_urine[i],"substrate");
					for (String key : single_import_result.keySet()) {
						if(single_import_result.get(key) == "substrate") {
							System.out.println(String.format("%s exerted to urine by %s.", enter_kidney_tubules_result.get(i).SMILES, enter_urine[i]));
							export_to_urine_flag = 1;
						}
				    }
				}
				
				if(export_to_urine_flag == 1) {
					System.out.println(String.format("%s exerted to urine.", enter_kidney_tubules_result.get(i).SMILES));
				}else {
					System.out.println(String.format("%s couldn't exert to urine.", enter_kidney_tubules_result.get(i).SMILES));
				}
			}
		}
		
		return UrineMetabolite;
		
		
		
		
	}
	
	/**
	 * Input: all possible generated metabolites: some produced metabolite will inhibit the cyp enzyme
	 * @param input
	 * @throws Exception 
	 */
	public static void InhibitionTestPhaseI(ArrayList<MetaboliteObject> input) throws Exception {
		// Inhibition testing I 
		// Test for all metabolite
		String[] cyp_enzyme_list = new String[] {"CYP1A2","CYP2B6","CYP2A6","CYP2C8","CYP2C9","CYP2C19","CYP2D6","CYP2E1","CYP3A4"};
		for(int i = 0; i< input.size(); i++) {
			Instances single_inhibitor_instance = RunClassification.generate_test_instance(input.get(i).SMILES, "inhibitor");
			for(int m = 0; m <cyp_enzyme_list.length; m++) {
				HashMap<String,String> single_inhibitor_result = phase1inhibitor.run_classifier(single_inhibitor_instance,cyp_enzyme_list[m]);
				for (String key : single_inhibitor_result.keySet()) {
					if(single_inhibitor_result.get(key) == "inhibitor") {
						input.get(i).AddTargetAsInhibitor(cyp_enzyme_list[m]);
						System.out.println(String.format("%s inhibit %s.", input.get(i).SMILES, cyp_enzyme_list[m]));
					}
			    }
			}
			
		}
	}
	
	
	/**
	 * input is all the generated metabolites + input molecules
	 * provide information for other phase II enzyme
	 * @param input
	 */
	public static void InhibitionTestPhaseII(ArrayList<MetaboliteObject> input) {
		
		
		
		
		
	}
	
	
	/**
	 * input is all the generated metabolites + input molecules
	 * provide information for other transporter on other location
	 * @param input
	 */
	public static void InhibitionTestTransporter(ArrayList<MetaboliteObject> input) {
		
		
		
		
		
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
