package xuan.drug_porter.som_phaseI;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Properties;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.io.listener.PropertiesListener;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import weka.classifiers.Classifier;
import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import xuan.drug_porter.descriptorUtils.GetAtomicDescriptors;
import xuan.drug_porter.descriptorUtils.SdfToSample;


public class SomPrediction {

	
	private static String current_dir = System.getProperty("user.dir");
	private static HashMap<String, String> som_model_path = CypModelPath();
	
	/**
	 * generate the test instance for predicting som
	 * @param input
	 * @return
	 * @throws CDKException 
	 * @throws FileNotFoundException 
	 */
	public static Instances create_test_instance(String input) throws CDKException, IOException {
		
		
		int nearest_atom = 3;
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer original_mole = builder.newInstance(IAtomContainer.class);
//		HashMap<IAtomContainer, Instances> hash_instance = new HashMap<IAtomContainer, Instances>();
//		ArrayList<Instance> instance_list = new ArrayList<Instance>();
		
	    // check if the input is sdf // or smiles;
	    if(input.contains(".sdf") || input.contains(".mol")) {
	    		original_mole = read_SDF_file(input);
	    }else {
	    		SmilesParser temp_smiles = new SmilesParser(builder);
		 	IAtomContainer atom_container = temp_smiles.parseSmiles(input);
		 	AtomContainerManipulator.suppressHydrogens(atom_container);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(atom_container);
		
		 	StructureDiagramGenerator sdg = new StructureDiagramGenerator();
			sdg.setMolecule(atom_container);
			sdg.generateCoordinates();
			original_mole = sdg.getMolecule();
	    }
//	    System.out.println(mole.getAtomCount()); //read smiles string
	    // add hydrogen 
//	    AtomContainerManipulator.convertImplicitToExplicitHydrogens(mole);
	    IAtomContainer mole = AtomContainerManipulator.removeHydrogens(original_mole);
	    System.out.println(mole.getAtomCount());
	    ArrayList<Attribute> attribute_name = generate_attribute_name(nearest_atom+1, 29); 
	     
	    FastVector<String> association = new FastVector<String>();
	    association.addElement("Yes");
		association.addElement("No");
		Attribute class_attribute = new Attribute("Class",association);
		attribute_name.add(class_attribute);
		 
		Instances test_instance = new Instances("Rel",attribute_name,mole.getAtomCount());
		
		test_instance.setClassIndex(class_attribute.index());
	    
		// getNearestAtoms may only works for 3d mol structure
	    ArrayList<ArrayList<String>> all_nearest_atom_set = GetAtomicDescriptors.getNearestAtoms(mole);
	    
	    int num_atom = mole.getAtomCount();
		for(int atoms = 0; atoms < num_atom; atoms++) {
			
			
		}
		
		for(int k = 0; k< num_atom; k++) {

			
			// iterate each instance;
			// convert each instance to weka instance;
			
			ArrayList<String> nearest_atom_set = all_nearest_atom_set.get(k);   // get index of atom but need to minus 1 because of chemsketch
			List<IAtom> atoms_set = new ArrayList<IAtom>();
			
			// add the # nearest atom into atom list for extract the atom descriptor;
			atoms_set.add(mole.getAtom(k));  									// add the original atoms
			for(int nna = 0; nna < nearest_atom; nna++) {
				IAtom tmp_atom = mole.getAtom(Integer.parseInt(nearest_atom_set.get(nna)));
				atoms_set.add(tmp_atom);
			}
			
			ArrayList<Double[]> descriptor_value = GetAtomicDescriptors.getAtomicDescriptor(mole,atoms_set, "");

			
			ArrayList<String> single_instance_value = new ArrayList<String>();
			for(int dv = 0; dv < descriptor_value.size(); dv++) {
				for (int dv2 = 0; dv2<descriptor_value.get(dv).length; dv2++) {
					single_instance_value.add(Double.toString(descriptor_value.get(dv)[dv2]));
				}
			}
			
			single_instance_value.add("No");
			
			
			// feature_set is for each temp feature.
			int num_of_attribute = test_instance.numAttributes();
			Instance feature_set = new DenseInstance(num_of_attribute); 
    		 
    		 	for(int f=0; f < num_of_attribute; f++) {
    		 		Attribute tmp_attr = attribute_name.get(f);
    		 		if(tmp_attr.isNumeric()) {
    		 			feature_set.setValue(tmp_attr,Double.parseDouble(single_instance_value.get(f)));
//    		 			feature_set.setValue(tmp_attr, single_instance_value.get(f));
    		 		}else if (tmp_attr.isNominal()) {
    		 			feature_set.setValue(tmp_attr, single_instance_value.get(f));
    		 		}
    			 
    		 	}
    		 
    		 	test_instance.add(feature_set);
			
		}
		test_instance.instance(0).setClassMissing();
//		System.out.println(test_instance.instance(0).toString());
//		System.exit(0);
		return test_instance;
		
		
	}
	
	
	/**
	 * 
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
	
	/**
	 * 
	 * @param factor
	 * @param num_of_descriptor
	 * @return
	 */
	public static ArrayList<Attribute> generate_attribute_name(int factor, int num_of_descriptor){
		
		ArrayList<Attribute> attribute = new ArrayList<Attribute>();
		int total_attribute_number = factor * num_of_descriptor;
		
		for(int i = 0; i < total_attribute_number; i++) {
			String each_attri = String.format("Attribute_%d", i);
			Attribute tmp_attribute = new Attribute(each_attri);
			attribute.add(tmp_attribute);
		}
		
		return attribute;
		
	}
	
	
	/**
	 * run prediction to get the predicted som
	 * model_type can only be string of "CYP1A2","CYP2A6","CYP2B6","CYP2C8","CYP2C9","CYP2C19","CYP2D6","CYP2E1","CYP3A4"
	 * should return list of site for appending new functional group
	 * @param test_instance
	 * @param model_type
	 * @throws Exception 
	 */
	public static HashMap<Integer,String> runSomClassifier(Instances instance, String model_type) throws Exception {
		
		HashMap<Integer,String> result = new HashMap<Integer,String>();
		
		if(model_type == "CYP1A2") {
			Classifier model = (Classifier) weka.core.SerializationHelper.read(som_model_path.get(model_type));
			result = runSomClassifierHelper(model,instance);
		}
		else if (model_type == "CYP2A6") {
			Classifier model = (Classifier) weka.core.SerializationHelper.read(som_model_path.get(model_type));
			result = runSomClassifierHelper(model,instance);
		}
		else if (model_type == "CYP2B6") {
			Classifier model = (Classifier) weka.core.SerializationHelper.read(som_model_path.get(model_type));
			result = runSomClassifierHelper(model,instance);
		}
		else if (model_type == "CYP2C8") {
			Classifier model = (Classifier) weka.core.SerializationHelper.read(som_model_path.get(model_type));
			result = runSomClassifierHelper(model,instance);
		}
		else if (model_type == "CYP2C9") {
			Classifier model = (Classifier) weka.core.SerializationHelper.read(som_model_path.get(model_type));
			result = runSomClassifierHelper(model,instance);
		}
		else if (model_type == "CYP2C19") {
			Classifier model = (Classifier) weka.core.SerializationHelper.read(som_model_path.get(model_type));
			result = runSomClassifierHelper(model,instance);
		}
		else if (model_type == "CYP2D6") {
			Classifier model = (Classifier) weka.core.SerializationHelper.read(som_model_path.get(model_type));
			result = runSomClassifierHelper(model,instance);
		}
		else if (model_type == "CYP2E1") {
			Classifier model = (Classifier) weka.core.SerializationHelper.read(som_model_path.get(model_type));
			result = runSomClassifierHelper(model,instance);
		}
		else if (model_type == "CYP3A4") {
			Classifier model = (Classifier) weka.core.SerializationHelper.read(som_model_path.get(model_type));
			result = runSomClassifierHelper(model,instance);
		}
		else {
			System.out.println("Unknown model type. (Please select from {\"CYP1A2\",\"CYP2A6\",\"CYP2B6\",\"CYP2C8\",\"CYP2C9\",\"CYP2C19\",\"CYP2D6\",\"CYP2E1\",\"CYP3A4\"}");
			System.exit(-1);
		}
		
		return result;
		
		
		
	}
	
	public static HashMap<Integer, String> runSomClassifierHelper(Classifier model, Instances data_instance) throws Exception {
		
		HashMap<Integer,String> result = new HashMap<Integer,String>();
		
		for(int i=0; i< data_instance.size(); i++) {
			double single_result = model.classifyInstance(data_instance.instance(i));
			if(single_result == 1.0) {
				result.put(i, "No");
			}else {
				result.put(i, "Yes");
			}
		}
		return result;
		
	}
	
	/**
	 * simple for building the path map
	 * @return
	 */
	public static HashMap<String, String> CypModelPath(){
		HashMap<String, String> maps = new HashMap<String,String>();
		maps.put("CYP1A2", String.format("%s/Model/PhaseISOMModel/%s_RANDOMFOREST.model", current_dir,"CYP1A2"));
		maps.put("CYP2A6", String.format("%s/Model/PhaseISOMModel/%s_RANDOMFOREST.model", current_dir,"CYP2A6"));
		maps.put("CYP2B6", String.format("%s/Model/PhaseISOMModel/%s_RANDOMFOREST.model", current_dir,"CYP2B6"));
		maps.put("CYP2C8", String.format("%s/Model/PhaseISOMModel/%s_RANDOMFOREST.model", current_dir,"CYP2C8"));
		maps.put("CYP2C9", String.format("%s/Model/PhaseISOMModel/%s_RANDOMFOREST.model", current_dir,"CYP2C9"));
		maps.put("CYP2C19", String.format("%s/Model/PhaseISOMModel/%s_RANDOMFOREST.model", current_dir,"CYP2C19"));
		maps.put("CYP2D6", String.format("%s/Model/PhaseISOMModel/%s_RANDOMFOREST.model", current_dir,"CYP2D6"));
		maps.put("CYP2E1", String.format("%s/Model/PhaseISOMModel/%s_RANDOMFOREST.model", current_dir,"CYP2E1"));
		maps.put("CYP3A4", String.format("%s/Model/PhaseISOMModel/%s_RANDOMFOREST.model", current_dir,"CYP3A4"));
		
		return maps;
		
	}
	
	
	public static void main(String[] args) throws Exception {
		String test_smiles = "/Users/xuan/Desktop/HMDB00001.sdf";
		Instances test_instance =  create_test_instance(test_smiles);
		HashMap<Integer,String> result = runSomClassifier(test_instance,"CYP1A2");
		for (Integer key : result.keySet()) {
			String result_each = result.get(key).toString();
			System.out.println(Integer.toString(key) +":"+ result_each);
	    }
		
		
	}
}
