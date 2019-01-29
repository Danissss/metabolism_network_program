package xuan.drug_porter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Properties;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.io.listener.PropertiesListener;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import com.opencsv.CSVReader;
import com.opencsv.CSVWriter;

import weka.classifiers.Classifier;
import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.ArffSaver;
import weka.core.converters.CSVLoader;
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
	
	
	
	protected static String current_dir = System.getProperty("user.dir");
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
	
	
	public static LinkedHashMap<String,String> run_classifier(String input, String transporter_name) throws Exception{
		
		LinkedHashMap<String,String> classified_result = new LinkedHashMap<String,String>();
		// list of model path; should around ~20 * 2 (both inhibitor and substrate model)
		String MDR1_model_inhibitor = current_dir +"/Model/MDR1_RANDOMFOREST_INHIBITOR.model";
		String MDR1_model_substrate = current_dir +"/Model/MDR1_RANDOMFOREST_SUBSTRATE.model";
		
		Instances substrate_instance = generate_test_instance(input,"substrate");
		Instances inhibitor_instance = generate_test_instance(input,"inhibitor");
//		System.out.print(transporter_name);
		if(transporter_name.contains("MDR1")) {
//			System.out.print(transporter_name);
			Classifier MDR1_substrate_model = (Classifier) weka.core.SerializationHelper.read(MDR1_model_substrate);
			double substrate_result = MDR1_substrate_model.classifyInstance(substrate_instance.get(0));
//			System.out.println(Double.toString(substrate_result));
			Classifier MDR1_inhibitor_model = (Classifier) weka.core.SerializationHelper.read(MDR1_model_inhibitor);
			double inhibitor_result = MDR1_inhibitor_model.classifyInstance(inhibitor_instance.get(0));
			
//			System.out.println(Double.toString(inhibitor_result));
			if(substrate_result == 0.0) {
				classified_result.put("MDR1_substrate", "non-substrate");
			}
			else {
				classified_result.put("MDR1_substrate", "substrate");
			}
			if(inhibitor_result == 0.0) {
				classified_result.put("MDR1_inhibitor", "non-inhibitor");
			}else {
				classified_result.put("MDR1_inhibitor", "inhibitor");
			}
			
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
	
	/*
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
			LinkedHashMap<String,String> result = run_classifier(input,transporter_name);
			for (String key : result.keySet()) {
				String result_each = result.get(key).toString();
				String result_key = key;
				System.out.println(result_key +":"+ result_each);
		    }
		}
		
		
		
	}
}
