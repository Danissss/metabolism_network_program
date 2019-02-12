package xuan.drug_porter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Properties;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.io.listener.PropertiesListener;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import com.google.common.collect.Maps;

import weka.classifiers.Classifier;
import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import xuan.drug_porter.descriptorUtils.FeatureGeneration;
import xuan.drug_porter.descriptorUtils.GenerateFeatureSingle;
import xuan.drug_porter.descriptorUtils.GeneratingFeatures;

/*
 * This function should take model and run classifiation for all transporters
 * 
 * 
 * 
 */
public class RunClassification {
	
	
	protected static String[] current_gene = new String[] {"ABCB11","ABCB1","ABCC2","ABCC3","ABCC4","ABCC5","ABCC6",
			"ABCG2","SLC10A2","SLC15A1","SLC15A2","SLC22A11","SLC22A12","SLC22A1","SLC22A2","SLC22A4",
			"SLC22A5","SLC22A6","SLC22A7","SLC47A1","SLC47A2","SLC51A","SLCO1A2","SLCO1B1",
			"SLCO1B3","SLCO2A1","SLCO2B1","SLCO4C1"};
	
	protected static String[] current_protein = new String[] {"BSEP","MDR1","MRP2","MRP3","MRP4","MRP5","MRP6",
			"BCRP1","ASBT","PEPT1","PEPT2","OAT4","URAT1","OCT1","OCT2","OCTN1",
			"OCTN2","OAT1","OAT2","MATE1","MATE2","OSTA","OATP1","OATP1B1",
			"OATP1B3","OATP2A1","OATP2B1","OATP4C1"};
	
	
	protected static String current_dir = System.getProperty("user.dir");
	private static HashMap<String, String> substrate_model_path = create_transporter_model_path("substrate");
	private static HashMap<String, String> inhibitor_model_path = create_transporter_model_path("inhibitor");
	private String sep = File.pathSeparator;
	/*
	 * ABCG2=BCRP=BCRP1: export: intestine_to_blood, brain_blood_barrier, liver_to_bile
	 * OATP=SLCO1A2: import: intestine_to_blood,
	 */
	/*
	 * transporter_name: by given the transporter name; return the transporter classified result
	 * if transporter_name is ALL, then calculate everything and return the result
	 * model_type = {inhibitor {inhibitor vs non-inhibitor}}{substrate {substrate vs non-substrate}}
	 */
	public static Instances generate_test_instance(String smiles, String model_type) throws Exception {
			 	
		FeatureGeneration featureGeneration = new FeatureGeneration();
	    GenerateFeatureSingle GFS = new GenerateFeatureSingle();
	    GeneratingFeatures GF = new GeneratingFeatures();
	    IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
	    IAtomContainer mole = builder.newInstance(IAtomContainer.class);
	    
	    // check if the input is sdf // or smiles;
	    if(smiles.contains(".sdf") || smiles.contains(".mol")) {
	    		mole = read_SDF_file(smiles);
	    }else {
	    		SmilesParser temp_smiles = new SmilesParser(builder);
		 	IAtomContainer atom_container = temp_smiles.parseSmiles(smiles);
		 	AtomContainerManipulator.suppressHydrogens(atom_container);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(atom_container);
		
		 	StructureDiagramGenerator sdg = new StructureDiagramGenerator();
			sdg.setMolecule(atom_container);
			sdg.generateCoordinates();
			mole = sdg.getMolecule();
	    }
	    
	 	
		 CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mole.getBuilder());
		
		 ArrayList<Attribute> attribute_al = GF.generateAllAttributes(mole);
	     
	     FastVector<String> association = new FastVector<String>();
		 if (model_type == "substrate") {
			association.addElement("substrate");
			association.addElement("non-substrate");
		 }
		 else {
			association.addElement("inhibitor");
			association.addElement("non-inhibitor");
		 }
		 Attribute class_attribute = new Attribute("Class",association);
		 attribute_al.add(class_attribute);
		 
		 Instances test_instance = new Instances("Rel",attribute_al,1);
		 test_instance.setClassIndex(class_attribute.index());

//		System.out.println(test_instance.numAttributes());
		String Features = GFS.generateOneinstance(mole,"ALL");
		int length = attribute_al.size(); //2282
		
		  String[] temp = Features.split(",");
//		  System.out.println(temp.length); //2282
		  Instance sample = new DenseInstance(length); 
		  for(int vidx = 0; vidx < temp.length; vidx++){
			  Attribute att = attribute_al.get(vidx);
			  
			  
			  String vle_string = temp[vidx];
			  if(vle_string.isEmpty()) {
				  System.out.println(att);
				  System.out.println(vle_string);
			  }

			  double vle = Double.parseDouble(temp[vidx]);
			  sample.setValue(att, vle);			

		  }
		  // it doesn't matter if set the class_attribute here; gonna change later
		  if(model_type == "substrate") {
			  sample.setValue(class_attribute, "non-substrate");
		  }else {
			  sample.setValue(class_attribute, "non-inhibitor");
		  }
		  
		  
		  test_instance.add(sample);
	    


		return test_instance;
	}
	
	
	/**
	 * return all transporter model result
	 * @param instance
	 * @return
	 * @throws Exception 
	 */
	
	public static HashMap<String,String> run_classifier_for_all_transporter(Instances instance) throws Exception {
		
		HashMap<String,String> result =  new HashMap<String,String>();
		String[] model_name = new String[] {"ABCB1"};
		
		for(int i=0 ; i<model_name.length; i++) {
			String s_model_name = model_name[i];
			HashMap<String,String> getback_result = run_classifier(instance,s_model_name);
			for (Map.Entry<String,String> entry : getback_result.entrySet()){
				result.put(entry.getKey(), entry.getValue());
			}
		}
		
		
		return result;
		
		
	}
	
	/**
	 * Get single transporter model 
	 * @param input
	 * @param transporter_name
	 * @return
	 * @throws Exception
	 */
	public static HashMap<String,String> run_classifier(Instances instance, String transporter_name) throws Exception{
		
		HashMap<String,String> classified_result = new HashMap<String,String>();
		
		Classifier MDR1_substrate_model = (Classifier) weka.core.SerializationHelper.read(substrate_model_path.get(transporter_name));
		double substrate_result = MDR1_substrate_model.classifyInstance(instance.get(0));
		Classifier MDR1_inhibitor_model = (Classifier) weka.core.SerializationHelper.read(inhibitor_model_path.get(transporter_name));
		double inhibitor_result = MDR1_inhibitor_model.classifyInstance(instance.get(0));
		String substrate = String.format("%s_substrate", transporter_name);
		String inhibitor = String.format("%s_inhibitor", transporter_name);
		if(substrate_result == 0.0) {
			classified_result.put(substrate, "non-substrate");
		}
		else {
			classified_result.put(substrate, "substrate");
		}
		if(inhibitor_result == 0.0) {
			classified_result.put(inhibitor, "non-inhibitor");
		}else {
			classified_result.put(inhibitor, "inhibitor");
		}
		
		return classified_result;
		
	}
	
	
	
	/**
	 * Read the sdf file and generate a IAtomContainerSet that contains all molecules in it.
	 * @param pathToInputFile
	 * @return
	 * @throws FileNotFoundException
	 * @throws CDKException
	 */

	public static IAtomContainer read_SDF_file(String pathToInputFile)
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
		
		IAtomContainer mole = MOLS.getAtomContainer(0);
		return mole;

	}
	
	
	public static HashMap<String,String> create_transporter_model_path(String sub_or_inh){
		HashMap<String,String> model_map = new HashMap<String,String>();
		if (sub_or_inh == "substrate") {
			for(int i = 0; i < current_gene.length; i++) {
				model_map.put(current_gene[i],String.format("%s/Model/TransporterModel/%s_RANDOMFOREST_SUBSTRATE.model", current_dir,current_gene[i]));
			}
		}
		else if(sub_or_inh == "inhibitor") {
			for(int i = 0; i < current_gene.length; i++) {
				model_map.put(current_gene[i],String.format("%s/Model/TransporterModel/%s_RANDOMFOREST_INHIBITOR.model", current_dir,current_gene[i]));
			}
		}
		
		return model_map;

	}
	
	/**
	 * create dictionary for gene name to map to protein name 
	 * @return
	 */
	public static HashMap<String,String> gene_name_protein_name_map(){
		
		HashMap<String, String> maps = new HashMap<String, String>();
		for(int i = 0; i < current_gene.length; i++) {
			maps.put(current_gene[i], current_protein[i]);
		}
		return maps;
	}
	
	
//	 {"ABCB11","ABCB1","ABCC2","ABCC3","ABCC4","ABCC5","ABCC6",
//			"ABCG2","SLC10A2","SLC15A1","SLC15A2","SLC22A11","SLC22A12","SLC22A1","SLC22A2","SLC22A4",
//			"SLC22A5","SLC22A6","SLC22A7","SLC47A1","SLC47A2","SLC51A","SLCO1A2","SLCO1B1",
//			"SLCO1B3","SLCO2A1","SLCO2B1","SLCO4C1"};
	
	/**
	 *  BBB = Brain-blood barrier
	 *  Naming convention: celltype_role_where
	 *  for example: Intestinal_Export_Intestine means on intestinal cell, export the compound, to intestine (i.e. export from blood to intestine)
	 * @return
	 */
	public static HashMap<String, String[]> protein_location(){
		HashMap<String, String[]> location_map = new HashMap<String, String[]>();
		location_map.put("ABCB11", new String[]{"Intestinal_Export_Intestine", "Hepatocytes_Export_Bile", "Kidney_Export_Urine", "BBB_Export_Blood"});
		location_map.put("ABCB1", new String[] {"Intestinal_Export_Intestine","Hepatocytes_Export_Bile","Kidney_Export_Urine","BBB_Export_Blood"});
		location_map.put("ABCC2", new String[] {"Intestinal_Export_Intestine","Hepatocytes_Export_Bile","Kidney_Export_Urine"});
		location_map.put("ABCC3", new String[] {""});
		location_map.put("ABCC4", new String[] {""});
		location_map.put("ABCC5", new String[] {""});
		location_map.put("ABCC6", new String[] {""});
		location_map.put("ABCG2", new String[] {""});
		location_map.put("SLC10A2", new String[] {""});
		location_map.put("SLC15A1", new String[] {""});
		location_map.put("SLC15A2", new String[] {""});
		location_map.put("SLC22A11", new String[] {""});
		location_map.put("SLC22A12", new String[] {""});
		location_map.put("SLC22A1", new String[] {""});
		location_map.put("SLC22A2", new String[] {""});
		location_map.put("SLC22A4", new String[] {""});
		location_map.put("SLC22A5", new String[] {""});
		location_map.put("SLC22A6", new String[] {""});
		location_map.put("SLC22A7", new String[] {""});
		location_map.put("SLC47A1", new String[] {""});
		location_map.put("SLC47A2", new String[] {""});
		location_map.put("SLC51A", new String[] {""});
		location_map.put("SLCO1A2", new String[] {""});
		location_map.put("SLCO1B1", new String[] {""});
		location_map.put("SLCO1B3", new String[] {""});
		location_map.put("SLCO2A1", new String[] {""});
		location_map.put("SLCO2B1", new String[] {""});
		location_map.put("SLCO4C1", new String[] {""});
		
		
		return location_map;
	}
	
	
	
	/**
	 * argument plan: either single transporter (by providing names) or all of them (by poviding "ALL" as arugment)
	 * smiles
	 * 
	 */
	public static void main(String[] args ) throws Exception {
		
		if(args.length < 2) {
			System.out.println("Too little arguments");
			System.exit(0);
		}
		
		// validation of smiles
		
		if(args[1] == "ALL") {
			// predict all transporter
			
		}
		else {
			String transporter_name = args[1];
			String input = args[0];
			Instances instance = generate_test_instance(input,"substrate");
			HashMap<String,String> result = run_classifier(instance,transporter_name);
			for (String key : result.keySet()) {
				String result_each = result.get(key).toString();
				String result_key = key;
				System.out.println(result_key +":"+ result_each);
		    }
		}
		
		
		
	}
}
