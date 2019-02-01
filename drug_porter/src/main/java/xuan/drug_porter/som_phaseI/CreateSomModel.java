package xuan.drug_porter.som_phaseI;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Random;

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
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import weka.attributeSelection.Ranker;
import weka.classifiers.evaluation.Evaluation;
import weka.classifiers.trees.RandomForest;
import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.ArffSaver;
import xuan.drug_porter.descriptorUtils.SdfToSample;
import xuan.drug_porter.descriptorUtils.GetAtomicDescriptors;

/**
 * 
 * create som predictor model for phaseI and phaseII
 * @author xuan
 *
 */
public class CreateSomModel {

	

	private static String current_dir = System.getProperty("user.dir");
	private static String[] enzyme_name = new String[] {"1A2","2A6","2B6","2C8","2C9","2C19","2D6","2E1","3A4"}; 
	
	public static ArrayList<Instances> create_phaseI_som_instance() throws CDKException, IOException {
		
		ArrayList<Instances> all_som_instance = new ArrayList<Instances>();
		ArrayList<ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>>> dataset =  create_phaseI_som_arraylist();
		int nearest_atom = 3;	// find nearest 4 atoms
		// previous step could be done by just call create_phaseI_som3d_arraylist() function
		// create instance;
		
		
		
		
	    
		 
		for(int i = 0; i< dataset.size(); i++) {
			String arff_path = current_dir+"/tmp_file/"+enzyme_name[i]+"_som.arff";
			
			
			// create the instance and attribute
			ArrayList<Attribute> attribute_name = generate_attribute_name(nearest_atom+1, 29); // return 4*29=116 attribute name;
			FastVector<String> association = new FastVector<String>();
		    association.addElement("Yes");
			association.addElement("No");
			Attribute class_attribute = new Attribute("Class",association);
			attribute_name.add(class_attribute);
			Instances new_instance = new Instances("Rel",attribute_name,10000);
			new_instance.setClassIndex(class_attribute.index()); 
			// create the instance and attribute
			int num_of_attribute =  new_instance.numAttributes();
			
			
			// iterate each enzyme 
			ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>> som_array = dataset.get(i);
			for(int j = 0; j< som_array.size(); j++) {
				
				
				
				// iterate each associated compounds
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> row = som_array.get(j);
				for (Map.Entry<IAtomContainer, ArrayList<ArrayList<String>>> entry : row.entrySet()){
					IAtomContainer mole = entry.getKey();
//					System.out.println(mole.getTitle());
					int num_atom = mole.getAtomCount();
					ArrayList<ArrayList<String>> all_nearest_atom_set = GetAtomicDescriptors.getNearestAtoms(mole);  // return nearest atom list in atom order
					ArrayList<ArrayList<String>> instances = entry.getValue();
//					instances contain list of instance; each instance contain atom index and yes or no
//					System.out.println(instances.size());
//					for each instance, 
					for(int k = 0; k<instances.size(); k++) {

						
						// iterate each instance;
						// convert each instance to weka instance;
						int som = Integer.parseInt(instances.get(k).get(0));
						String is_som = instances.get(k).get(1);
						ArrayList<String> nearest_atom_set = all_nearest_atom_set.get(som-1);   // get index of atom but need to minus 1 because of chemsketch
						List<IAtom> atoms_set = new ArrayList<IAtom>();
						
						// add the # nearest atom into atom list for extract the atom descriptor;
						atoms_set.add(mole.getAtom(som-1));  									// add the original atoms
						for(int nna = 0; nna < nearest_atom; nna++) {
							IAtom tmp_atom = mole.getAtom(Integer.parseInt(nearest_atom_set.get(nna)));
							atoms_set.add(tmp_atom);
						}
						
						ArrayList<Double[]> descriptor_value = GetAtomicDescriptors.getAtomicDescriptor(mole,atoms_set, "");
						
						// concatate the double[];
						// add the true value;
						// make it weka instance;
						ArrayList<String> single_instance_value = new ArrayList<String>();
						for(int dv = 0; dv < descriptor_value.size(); dv++) {
							for (int dv2 = 0; dv2<descriptor_value.get(dv).length; dv2++) {
								single_instance_value.add(Double.toString(descriptor_value.get(dv)[dv2]));
							}
						}
						
						single_instance_value.add(is_som);
						
						
						// feature_set is for each temp feature.
						Instance feature_set = new DenseInstance(num_of_attribute); // num_of_attribute=2281
			    		 
			    		 	for(int f=0; f < num_of_attribute; f++) {
			    		 		Attribute tmp_attr = attribute_name.get(f);
			    		 		if(tmp_attr.isNumeric()) {
			    		 			feature_set.setValue(tmp_attr,Double.parseDouble(single_instance_value.get(f)));
//			    		 			feature_set.setValue(tmp_attr, single_instance_value.get(f));
			    		 		}else if (tmp_attr.isNominal()) {
			    		 			feature_set.setValue(tmp_attr, single_instance_value.get(f));
			    		 		}
			    			 
			    		 	}
			    		 
			    		 	new_instance.add(feature_set);
						
					}
					
					
				}
				
			}
			
			
			File arff_instance = new File(arff_path);
			ArffSaver arffSaver = new ArffSaver();
			arffSaver.setInstances(new_instance);
			arffSaver.setFile(arff_instance);
			arffSaver.writeBatch();
			
			
			all_som_instance.add(new_instance);
			
			System.out.println(arff_path+" finished!======================!");
			
//			new_instance.clear();
//			attribute_name.clear();
		}
		
		return all_som_instance;
		
	}
	
	/**
	 * parse the smiles
	 * if the input is smiles string, parse it to mole
	 * if it is mol or sdf file, parse it to mole
	 * @param smiles
	 * @return
	 * @throws CDKException 
	 * @throws IOException 
	 */
	
	
	public static HashMap<IAtomContainer, Instances> create_test_instance(String input) throws CDKException, IOException {
		
		
		int nearest_atom = 3;
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer mole = builder.newInstance(IAtomContainer.class);
		HashMap<IAtomContainer, Instances> hash_instance = new HashMap<IAtomContainer, Instances>();
		ArrayList<Instance> instance_list = new ArrayList<Instance>();
		
	    // check if the input is sdf // or smiles;
	    if(input.contains(".sdf") || input.contains(".mol")) {
	    		mole = read_SDF_file(input);
	    }else {
	    		SmilesParser temp_smiles = new SmilesParser(builder);
		 	IAtomContainer atom_container = temp_smiles.parseSmiles(input);
		 	AtomContainerManipulator.suppressHydrogens(atom_container);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(atom_container);
		
		 	StructureDiagramGenerator sdg = new StructureDiagramGenerator();
			sdg.setMolecule(atom_container);
			sdg.generateCoordinates();
			mole = sdg.getMolecule();
	    }
	    
	    
	    ArrayList<Attribute> attribute_name = generate_attribute_name(nearest_atom+1, 29); 
	     
	    FastVector<String> association = new FastVector<String>();
	    association.addElement("Yes");
		association.addElement("No");
		Attribute class_attribute = new Attribute("Class",association);
		attribute_name.add(class_attribute);
		 
		Instances test_instance = new Instances("Rel",attribute_name,mole.getAtomCount());
		
		test_instance.setClassIndex(class_attribute.index());
	    
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
			
			single_instance_value.add("?");
			
			
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
		
		hash_instance.put(mole,test_instance);
	    
	    
	    
	    // return the HashMap<IAtomContainer, Instances>
		// each Instances contain descriptor value of each atom
		// in prediction
		// test if the instance has the site.
		
		return hash_instance;
		
		
	}
	
	
	
	
	
	/**
	 * why the clear() will clear the stuff even after insert into map
	 * dataList.add(map) will put a reference to map in the list, so it's not a copy. 
	 * When you then do map.clear() afterwards, it erases the content of the map in the list too, 
	 * because it is the very same object. Do dataList.add(map.clone()) instead or (preferably) 
	 * do map = new HashMap<>(); afterwards.
	 * 
	 * @return 3D arraylist
	 * @throws FileNotFoundException
	 * @throws CDKException
	 */
	public static ArrayList<ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>>> create_phaseI_som_arraylist() throws FileNotFoundException, CDKException {
		
		String som_sdf_file_path = current_dir+"/tmp_file/SOMDATA.sdf";
		ArrayList<ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>>> dataset = new ArrayList<ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>>>();
		
		//ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>> more like ArrayList<HashMap<IAtomContainer,ArrayList<String[]>>> 
		// ArrayList<String[]> store the number of instance from the IAtomContainer.
		ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>> hash_1A2 = new ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>>();
		ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>> hash_2A6 = new ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>>();
		ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>> hash_2B6 = new ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>>();
		ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>> hash_2C8 = new ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>>();
		ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>> hash_2C9 = new ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>>();
		ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>> hash_2C19 = new ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>>();
		ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>> hash_2D6 = new ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>>();
		ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>> hash_2E1 = new ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>>();
		ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>> hash_3A4 = new ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>>();
		

		// hashmap<mol,ArrayList<ArrayList<>>
		
		
		
		SdfToSample sdftosample = new SdfToSample();
		IAtomContainerSet mol_set = sdftosample.readFile(som_sdf_file_path);
		
		for(int i=0; i<mol_set.getAtomContainerCount();i++) {
			
			ArrayList<ArrayList<String>> list_1A2 = new ArrayList<ArrayList<String>>();
			ArrayList<ArrayList<String>> list_2A6 = new ArrayList<ArrayList<String>>();
			ArrayList<ArrayList<String>> list_2B6 = new ArrayList<ArrayList<String>>();
			ArrayList<ArrayList<String>> list_2C8 = new ArrayList<ArrayList<String>>();
			ArrayList<ArrayList<String>> list_2C9 = new ArrayList<ArrayList<String>>();
			ArrayList<ArrayList<String>> list_2C19 = new ArrayList<ArrayList<String>>();
			ArrayList<ArrayList<String>> list_2D6 = new ArrayList<ArrayList<String>>();
			ArrayList<ArrayList<String>> list_2E1 = new ArrayList<ArrayList<String>>();
			ArrayList<ArrayList<String>> list_3A4 = new ArrayList<ArrayList<String>>();
			
			IAtomContainer mol = mol_set.getAtomContainer(i);
			System.out.println(mol.getTitle());
//			CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mol.getBuilder());
			SmilesGenerator sg      = new SmilesGenerator(SmiFlavor.Isomeric);
			String smi_mol     = sg.create(mol);
			
			
			
			// this contain the right order to atoms 
			ArrayList<Integer> atom_list_1A2 = new ArrayList<Integer>();
			for(int atom_num = 0; atom_num < mol.getAtomCount(); atom_num++) {
				IAtom tmp_atom = mol.getAtom(atom_num);
				int num = mol.getAtomNumber(tmp_atom);
//				System.out.println(num+":"+tmp_atom.getSymbol());
				atom_list_1A2.add(num+1);
			}
			ArrayList<Integer> atom_list_2A6 = new ArrayList<Integer>(atom_list_1A2);
			ArrayList<Integer> atom_list_2B6 = new ArrayList<Integer>(atom_list_1A2);
			ArrayList<Integer> atom_list_2C8 = new ArrayList<Integer>(atom_list_1A2);
			ArrayList<Integer> atom_list_2C9 = new ArrayList<Integer>(atom_list_1A2);
			ArrayList<Integer> atom_list_2C19 = new ArrayList<Integer>(atom_list_1A2);
			ArrayList<Integer> atom_list_2D6 = new ArrayList<Integer>(atom_list_1A2);
			ArrayList<Integer> atom_list_2E1 = new ArrayList<Integer>(atom_list_1A2);
			ArrayList<Integer> atom_list_3A4 = new ArrayList<Integer>(atom_list_1A2);
			
			
			
			// This is for mark if the compound has the SOM; since, not all the compound is metabolized by all the enzyme
			// Can reduce the imbalanced dataset problem
			int flag_1A2 = 0;
			int flag_2A6 = 0;
			int flag_2B6 = 0;
			int flag_2C8 = 0;
			int flag_2C9 = 0;
			int flag_2C19 = 0;
			int flag_2D6 = 0;
			int flag_2E1 = 0;
			int flag_3A4 = 0;
			
			
			Map<Object,Object> prop =   mol.getProperties();
			for (Map.Entry<Object, Object> entry : prop.entrySet()){
//				System.out.println(entry.getKey() + "/" + entry.getValue());
				
				String key = entry.getKey().toString();
				
				
				
				
				
				
				
				
				
				
				
				if (key.contains("1A2")) {
					String values = entry.getValue().toString();
					String[] value_list = values.split(" ");
					for(int vl = 0; vl<value_list.length; vl++) {
						ArrayList<String> array_1A2 = new ArrayList<String>();
						int value_int = Integer.valueOf(value_list[vl]);
						atom_list_1A2.remove(Integer.valueOf(value_int)); // this line remove the element based on the value
						array_1A2.add(Integer.toString(value_int));
						array_1A2.add("Yes");
						list_1A2.add(array_1A2);		// add to the list_1A2 (ArrayList<ArrayList<String>>)
//						System.out.println(array_1A2.toString());
//						array_1A2.clear(); 			// clean the array_1A2; this will clear everything inside the list_1A2;
					}
					
					
					flag_1A2 = 1;
				}
				else if (key.contains("2A6")) {
					String values = entry.getValue().toString();
					String[] value_list = values.split(" ");
					for(int vl = 0; vl<value_list.length; vl++) {
						ArrayList<String> array_2A6 = new ArrayList<String>();
						int value_int = Integer.valueOf(value_list[vl]);
						atom_list_2A6.remove(Integer.valueOf(value_list[vl]));
						array_2A6.add(Integer.toString(value_int));
						array_2A6.add("Yes");
						list_2A6.add(array_2A6);
					}
					flag_2A6 = 1;
				}
				else if (key.contains("2B6")) {
					String values = entry.getValue().toString();
					String[] value_list = values.split(" ");
					for(int vl = 0; vl<value_list.length; vl++) {
						ArrayList<String> array_2B6 = new ArrayList<String>();
						int value_int = Integer.valueOf(value_list[vl]);
						atom_list_2B6.remove(Integer.valueOf(value_int));
						array_2B6.add(Integer.toString(value_int));
						array_2B6.add("Yes");
						list_2B6.add(array_2B6);
					}
					flag_2B6 = 1;
				}
				else if (key.contains("2C8")) {
					String values = entry.getValue().toString();
					String[] value_list = values.split(" ");
					for(int vl = 0; vl<value_list.length; vl++) {
						ArrayList<String> array_2C8 = new ArrayList<String>();
						int value_int = Integer.valueOf(value_list[vl]);
						atom_list_2C8.remove(Integer.valueOf(value_int));
						array_2C8.add(Integer.toString(value_int));
						array_2C8.add("Yes");
						list_2C8.add(array_2C8);
					}
					flag_2C8 = 1;
				}
				else if (key.contains("2C9")) {
					String values = entry.getValue().toString();
					String[] value_list = values.split(" ");
					for(int vl = 0; vl<value_list.length; vl++) {
						ArrayList<String> array_2C9 = new ArrayList<String>();
						int value_int = Integer.valueOf(value_list[vl]);
						atom_list_2C9.remove(Integer.valueOf(value_int));
						array_2C9.add(Integer.toString(value_int));
						array_2C9.add("Yes");
						list_2C9.add(array_2C9);
					}
					flag_2C9 = 1;
				}
				else if (key.contains("2C19")) {
					String values = entry.getValue().toString();
					String[] value_list = values.split(" ");
					for(int vl = 0; vl<value_list.length; vl++) {
						ArrayList<String> array_2C19 = new ArrayList<String>();
						int value_int = Integer.valueOf(value_list[vl]);
						atom_list_2C19.remove(Integer.valueOf(value_int));
						array_2C19.add(Integer.toString(value_int));
						array_2C19.add("Yes");
						list_2C19.add(array_2C19);
					}
					flag_2C19 = 1;
				}
				else if (key.contains("2D6")) {
					String values = entry.getValue().toString();
					String[] value_list = values.split(" ");
					for(int vl = 0; vl<value_list.length; vl++) {
						ArrayList<String> array_2D6 = new ArrayList<String>();
						int value_int = Integer.valueOf(value_list[vl]);
						atom_list_2D6.remove(Integer.valueOf(value_int));
						array_2D6.add(Integer.toString(value_int));
						array_2D6.add("Yes");
						list_2D6.add(array_2D6);
					}
					flag_2D6 = 1;
				}
				else if (key.contains("2E1")) {
					String values = entry.getValue().toString();
					String[] value_list = values.split(" ");
					for(int vl = 0; vl<value_list.length; vl++) {
						ArrayList<String> array_2E1 = new ArrayList<String>();
						int value_int = Integer.valueOf(value_list[vl]);
						atom_list_2E1.remove(Integer.valueOf(value_int));
						array_2E1.add(Integer.toString(value_int));
						array_2E1.add("Yes");
						list_2E1.add(array_2E1);
					}
					flag_2E1 = 1;
				}
				else if (key.contains("3A4")) {
					String values = entry.getValue().toString();
					String[] value_list = values.split(" ");
					for(int vl = 0; vl<value_list.length; vl++) {
						ArrayList<String> array_3A4 = new ArrayList<String>();
						int value_int = Integer.valueOf(value_list[vl]);
						atom_list_3A4.remove(Integer.valueOf(value_int));
						array_3A4.add(Integer.toString(value_int));
						array_3A4.add("Yes");
						list_3A4.add(array_3A4);
					}
					flag_3A4 = 1;
				}
				
				
			}
			
			// add negative instance.
			if(flag_1A2==1) {
				
				for (int left_atom = 0; left_atom < atom_list_1A2.size();left_atom++) {
					ArrayList<String> array_1A2 = new ArrayList<String>();
					array_1A2.add(Integer.toString(atom_list_1A2.get(left_atom)));
					array_1A2.add("No");
					list_1A2.add(array_1A2);
//					System.out.println(array_1A2.toString());
//					array_1A2.clear();
					
				}
				
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
//				System.out.println("list_1A2.size(): "+list_1A2.size());
				properties.put(mol,list_1A2);
				hash_1A2.add(properties);
				
//				for (Map.Entry<IAtomContainer, ArrayList<ArrayList<String>>> entry : properties.entrySet()){
//					System.out.println("hashmap");
//					ArrayList<ArrayList<String>> tmp = entry.getValue();
//					System.out.println(tmp.size());
//					
//					for(int mm = 0; mm < tmp.size(); mm++) {
//						System.out.println(tmp.get(mm).toString());
//					}
//					
//				}
//				System.out.println("====");
			}
			
			if(flag_2A6 == 1) {
				for (int left_atom = 0; left_atom < atom_list_2A6.size();left_atom++) {
					ArrayList<String> array_2A6 = new ArrayList<String>();
					array_2A6.add(Integer.toString(atom_list_2A6.get(left_atom)));
					array_2A6.add("No");
					list_2A6.add(array_2A6);
				}
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
				properties.put(mol,list_2A6);
				hash_2A6.add(properties);
			}
			
			if(flag_2B6 == 1) {
				for (int left_atom = 0; left_atom < atom_list_2B6.size();left_atom++) {
					ArrayList<String> array_2B6 = new ArrayList<String>();
					array_2B6.add(Integer.toString(atom_list_2B6.get(left_atom)));
					array_2B6.add("No");
					list_2B6.add(array_2B6);
				}
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
				properties.put(mol,list_2B6);
				hash_2B6.add(properties);
			}
			
			if(flag_2C8 == 1) {
				for (int left_atom = 0; left_atom < atom_list_2C8.size();left_atom++) {	
					ArrayList<String> array_2C8 = new ArrayList<String>();
					array_2C8.add(Integer.toString(atom_list_2C8.get(left_atom)));
					array_2C8.add("No");
					list_2C8.add(array_2C8);
				}
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
				properties.put(mol,list_2C8);
				hash_2C8.add(properties);
			}
			
			if(flag_2C9 == 1) {
				for (int left_atom = 0; left_atom < atom_list_2C9.size();left_atom++) {
					ArrayList<String> array_2C9 = new ArrayList<String>();
					array_2C9.add(Integer.toString(atom_list_2C9.get(left_atom)));
					array_2C9.add("No");
					list_2C9.add(array_2C9);
				}
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
				properties.put(mol,list_2C9);
				hash_2C9.add(properties);
			}
			
			if(flag_2C19 == 1) {
				for (int left_atom = 0; left_atom < atom_list_2C19.size();left_atom++) {
					ArrayList<String> array_2C19 = new ArrayList<String>();
					array_2C19.add(Integer.toString(atom_list_2C19.get(left_atom)));
					array_2C19.add("No");
					list_2C19.add(array_2C19);
				}
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
				properties.put(mol,list_2C19);
				hash_2C19.add(properties);
			}
			
			if(flag_2D6 == 1) {		
				for (int left_atom = 0; left_atom < atom_list_2D6.size();left_atom++) {
					ArrayList<String> array_2D6 = new ArrayList<String>();
					array_2D6.add(Integer.toString(atom_list_2D6.get(left_atom)));
					array_2D6.add("No");
					list_2D6.add(array_2D6);
				}
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
				properties.put(mol,list_2D6);
				hash_2D6.add(properties);
			}
			
			if(flag_2E1 == 1) {
				for (int left_atom = 0; left_atom < atom_list_2E1.size();left_atom++) {
					ArrayList<String> array_2E1 = new ArrayList<String>();
					array_2E1.add(Integer.toString(atom_list_2E1.get(left_atom)));
					array_2E1.add("No");
					list_2E1.add(array_2E1);
				}
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
				properties.put(mol,list_2E1);
				hash_2E1.add(properties);
			}
			
			if(flag_3A4 == 1) {	
				for (int left_atom = 0; left_atom < atom_list_3A4.size();left_atom++) {
					ArrayList<String> array_3A4 = new ArrayList<String>();
					array_3A4.add(Integer.toString(atom_list_3A4.get(left_atom)));
					array_3A4.add("No");
					list_3A4.add(array_3A4);
				}
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
				properties.put(mol,list_3A4);
				hash_3A4.add(properties);
			}
			
		}
		
//		System.out.println("hash_1A2.size(): "+hash_1A2.size());
//		HashMap<IAtomContainer,ArrayList<ArrayList<String>>> tmp_hash = hash_1A2.get(0);
//		for (Map.Entry<IAtomContainer, ArrayList<ArrayList<String>>> entry : tmp_hash.entrySet()){
//			System.out.println("hashmap");
//			ArrayList<ArrayList<String>> tmp = entry.getValue();
//			System.out.println(tmp.size());
//			
//			for(int mm = 0; mm < tmp.size(); mm++) {
//				System.out.println(tmp.get(mm).toString());
//			}
//			
//		}
//		for (Map.Entry<IAtomContainer,ArrayList<ArrayList<String>>> entry : hash_1A2.entrySet()){
//		
//		}
		
		dataset.add(hash_1A2);
		dataset.add(hash_2A6);
		dataset.add(hash_2B6);
		dataset.add(hash_2C8);
		dataset.add(hash_2C9);
		dataset.add(hash_2C19);
		dataset.add(hash_2D6);
		dataset.add(hash_2E1);
		dataset.add(hash_3A4);
		
		
		return dataset;
		
	}
	
	/**
	 * 
	 * 
	 * @return 3D arraylist
	 * @throws FileNotFoundException
	 * @throws CDKException
	 */
	public static ArrayList<ArrayList<ArrayList<String>>> create_phaseI_som_3d_arrayList() throws FileNotFoundException, CDKException {
		
		String som_sdf_file_path = current_dir+"/tmp_file/SOMDATA.sdf";
		ArrayList<ArrayList<ArrayList<String>>> dataset = new ArrayList<ArrayList<ArrayList<String>>>();
		
		ArrayList<ArrayList<String>> list_1A2 = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> list_2A6 = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> list_2B6 = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> list_2C8 = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> list_2C9 = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> list_2C19 = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> list_2D6 = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> list_2E1 = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> list_3A4 = new ArrayList<ArrayList<String>>();
		
		
		SdfToSample sdftosample = new SdfToSample();
		IAtomContainerSet mol_set = sdftosample.readFile(som_sdf_file_path);
		
		for(int i=0; i<mol_set.getAtomContainerCount();i++) {
			
			IAtomContainer mol = mol_set.getAtomContainer(i);
//			CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mol.getBuilder());
			SmilesGenerator sg      = new SmilesGenerator(SmiFlavor.Isomeric);
			String smi_mol     = sg.create(mol);
			
			
			ArrayList<Integer> atom_list_1A2 = new ArrayList<Integer>();
			for(int atom_num = 0; atom_num < mol.getAtomCount(); atom_num++) {
				IAtom tmp_atom = mol.getAtom(atom_num);
				int num = mol.getAtomNumber(tmp_atom);
//				System.out.println(num+":"+tmp_atom.getSymbol());
				atom_list_1A2.add(num+1);
			}
			ArrayList<Integer> atom_list_2A6 = new ArrayList<Integer>(atom_list_1A2);
			ArrayList<Integer> atom_list_2B6 = new ArrayList<Integer>(atom_list_1A2);
			ArrayList<Integer> atom_list_2C8 = new ArrayList<Integer>(atom_list_1A2);
			ArrayList<Integer> atom_list_2C9 = new ArrayList<Integer>(atom_list_1A2);
			ArrayList<Integer> atom_list_2C19 = new ArrayList<Integer>(atom_list_1A2);
			ArrayList<Integer> atom_list_2D6 = new ArrayList<Integer>(atom_list_1A2);
			ArrayList<Integer> atom_list_2E1 = new ArrayList<Integer>(atom_list_1A2);
			ArrayList<Integer> atom_list_3A4 = new ArrayList<Integer>(atom_list_1A2);
			
			
			
			// This is for mark if the compound has the SOM; since, not all the compound is metabolized by all the enzyme
			// Can reduce the imbalanced dataset problem
			int flag_1A2 = 0;
			int flag_2A6 = 0;
			int flag_2B6 = 0;
			int flag_2C8 = 0;
			int flag_2C9 = 0;
			int flag_2C19 = 0;
			int flag_2D6 = 0;
			int flag_2E1 = 0;
			int flag_3A4 = 0;
			
			
			Map<Object,Object> prop =   mol.getProperties();
			for (Map.Entry<Object, Object> entry : prop.entrySet()){
//				System.out.println(entry.getKey() + "/" + entry.getValue());
				
				String key = entry.getKey().toString();
				ArrayList<String> array_1A2 = new ArrayList<String>();
				ArrayList<String> array_2A6 = new ArrayList<String>();
				ArrayList<String> array_2B6 = new ArrayList<String>();
				ArrayList<String> array_2C8 = new ArrayList<String>();
				ArrayList<String> array_2C9 = new ArrayList<String>();
				ArrayList<String> array_2C19 = new ArrayList<String>();
				ArrayList<String> array_2D6 = new ArrayList<String>();
				ArrayList<String> array_2E1 = new ArrayList<String>();
				ArrayList<String> array_3A4 = new ArrayList<String>();
				
				
				if (key.contains("1A2")) {
					int value_int = Integer.valueOf(entry.getValue().toString());
					atom_list_1A2.remove(value_int-1);
					array_1A2.add(smi_mol);
					array_1A2.add(Integer.toString(value_int));
					array_1A2.add("Yes");
					list_1A2.add(array_1A2);
					array_1A2.clear();
					flag_1A2 = 1;
				}
				else if (key.contains("2A6")) {
					int value_int = Integer.valueOf(entry.getValue().toString());
					atom_list_2A6.remove(value_int-1);
					array_2A6.add(smi_mol);
					array_2A6.add(Integer.toString(value_int));
					array_2A6.add("Yes");
					list_2A6.add(array_2A6);
					array_2A6.clear();
					flag_2A6 = 1;
					
				}
				else if (key.contains("2B6")) {
					int value_int = Integer.valueOf(entry.getValue().toString());
					atom_list_2B6.remove(value_int-1);
					array_2B6.add(smi_mol);
					array_2B6.add(Integer.toString(value_int));
					array_2B6.add("Yes");
					list_2B6.add(array_2B6);
					array_2B6.clear();
					flag_2B6 = 1;
				}
				else if (key.contains("2C8")) {
					int value_int = Integer.valueOf(entry.getValue().toString());
					atom_list_2C8.remove(value_int-1);
					array_2C8.add(smi_mol);
					array_2C8.add(Integer.toString(value_int));
					array_2C8.add("Yes");
					list_2C8.add(array_2C8);
					array_2C8.clear();
					flag_2C8 = 1;
				}
				else if (key.contains("2C9")) {
					int value_int = Integer.valueOf(entry.getValue().toString());
					atom_list_2C9.remove(value_int-1);
					array_2C9.add(smi_mol);
					array_2C9.add(Integer.toString(value_int));
					array_2C9.add("Yes");
					list_2C9.add(array_2C9);
					array_2C9.clear();
					flag_2C9 = 1;
				}
				else if (key.contains("2C19")) {
					int value_int = Integer.valueOf(entry.getValue().toString());
					atom_list_2C19.remove(value_int-1);
					array_2C19.add(smi_mol);
					array_2C19.add(Integer.toString(value_int));
					array_2C19.add("Yes");
					list_2C19.add(array_2C19);
					array_2C19.clear();
					flag_2C19 = 1;
				}
				else if (key.contains("2D6")) {
					int value_int = Integer.valueOf(entry.getValue().toString());
					atom_list_2D6.remove(value_int-1);
					array_2D6.add(smi_mol);
					array_2D6.add(Integer.toString(value_int));
					array_2D6.add("Yes");
					list_2D6.add(array_2D6);
					array_2D6.clear();
					flag_2D6 = 1;
				}
				else if (key.contains("2E1")) {
					int value_int = Integer.valueOf(entry.getValue().toString());
					atom_list_2E1.remove(value_int-1);
					array_2E1.add(smi_mol);
					array_2E1.add(Integer.toString(value_int));
					array_2E1.add("Yes");
					list_2E1.add(array_2E1);
					array_2E1.clear();
					flag_2E1 = 1;
				}
				else if (key.contains("3A4")) {
					int value_int = Integer.valueOf(entry.getValue().toString());
					atom_list_3A4.remove(value_int-1);
					array_3A4.add(smi_mol);
					array_3A4.add(Integer.toString(value_int));
					array_3A4.add("Yes");
					list_3A4.add(array_3A4);
					array_3A4.clear();
					flag_3A4 = 1;
				}
				
				
			}
			
			// add negative instance.
			if(flag_1A2==1) {
				ArrayList<String> array_1A2 = new ArrayList<String>();
				for (int left_atom = 0; left_atom < atom_list_1A2.size();left_atom++) {
					array_1A2.add(smi_mol);
					array_1A2.add(Integer.toString(atom_list_1A2.get(left_atom)));
					array_1A2.add("No");
					list_1A2.add(array_1A2);
					array_1A2.clear();
				}
			}
			
			if(flag_2A6 == 1) {
				ArrayList<String> array_2A6 = new ArrayList<String>();
				for (int left_atom = 0; left_atom < atom_list_2A6.size();left_atom++) {
					array_2A6.add(smi_mol);
					array_2A6.add(Integer.toString(atom_list_2A6.get(left_atom)));
					array_2A6.add("No");
					list_2A6.add(array_2A6);
					array_2A6.clear();
				}
			}
			
			if(flag_2B6 == 1) {
				ArrayList<String> array_2B6 = new ArrayList<String>();
				for (int left_atom = 0; left_atom < atom_list_2B6.size();left_atom++) {
					array_2B6.add(smi_mol);
					array_2B6.add(Integer.toString(atom_list_2B6.get(left_atom)));
					array_2B6.add("No");
					list_2B6.add(array_2B6);
					array_2B6.clear();
				}
			}
			
			if(flag_2C8 == 1) {
				ArrayList<String> array_2C8 = new ArrayList<String>();
				for (int left_atom = 0; left_atom < atom_list_2C8.size();left_atom++) {
					array_2C8.add(smi_mol);
					array_2C8.add(Integer.toString(atom_list_2C8.get(left_atom)));
					array_2C8.add("No");
					list_2C8.add(array_2C8);
					array_2C8.clear();
				}
			}
			
			if(flag_2C9 == 1) {
				ArrayList<String> array_2C9 = new ArrayList<String>();
				for (int left_atom = 0; left_atom < atom_list_2C9.size();left_atom++) {
					array_2C9.add(smi_mol);
					array_2C9.add(Integer.toString(atom_list_2C9.get(left_atom)));
					array_2C9.add("No");
					list_2C9.add(array_2C9);
					array_2C9.clear();
				}
			}
			
			if(flag_2C19 == 1) {
				ArrayList<String> array_2C19 = new ArrayList<String>();
				for (int left_atom = 0; left_atom < atom_list_2C19.size();left_atom++) {
					array_2C19.add(smi_mol);
					array_2C19.add(Integer.toString(atom_list_2C19.get(left_atom)));
					array_2C19.add("No");
					list_2C19.add(array_2C19);
					array_2C19.clear();
				}
			}
			
			if(flag_2D6 == 1) {
				ArrayList<String> array_2D6 = new ArrayList<String>();
				for (int left_atom = 0; left_atom < atom_list_2D6.size();left_atom++) {
					array_2D6.add(smi_mol);
					array_2D6.add(Integer.toString(atom_list_2D6.get(left_atom)));
					array_2D6.add("No");
					list_2D6.add(array_2D6);
					array_2D6.clear();
				}
			}
			
			if(flag_2E1 == 1) {
				ArrayList<String> array_2E1 = new ArrayList<String>();
				for (int left_atom = 0; left_atom < atom_list_2E1.size();left_atom++) {
					array_2E1.add(smi_mol);
					array_2E1.add(Integer.toString(atom_list_2E1.get(left_atom)));
					array_2E1.add("No");
					list_2E1.add(array_2E1);
					array_2E1.clear();
				}
			}
			
			if(flag_3A4 == 1) {
				ArrayList<String> array_3A4 = new ArrayList<String>();
				for (int left_atom = 0; left_atom < atom_list_3A4.size();left_atom++) {
					array_3A4.add(smi_mol);
					array_3A4.add(Integer.toString(atom_list_3A4.get(left_atom)));
					array_3A4.add("No");
					list_3A4.add(array_3A4);
					array_3A4.clear();
				}
			}
		}	
		
		dataset.add(list_1A2);
		dataset.add(list_2A6);
		dataset.add(list_2B6);
		dataset.add(list_2C8);
		dataset.add(list_2C9);
		dataset.add(list_2C19);
		dataset.add(list_2D6);
		dataset.add(list_2E1);
		dataset.add(list_3A4);
		
		
		return dataset;
		
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
	 * two way to create model
	 * run the new SOM sdf file and generate new instance
	 * or load the existing arff file 
	 * @param all_som_instances
	 * @throws Exception 
	 */
	public static void create_som_model(ArrayList<Instances> all_som_instances, String option) throws Exception {
		
		
		for(int s = 0; s < all_som_instances.size(); s++ ) {
			String model_path =  current_dir+"/tmp_file/"+enzyme_name[s]+"RANDONFOREST_ATOMIC.model";
			Instances single_instance = all_som_instances.get(s);
			
			
			int num_attribute = single_instance.numAttributes();
			for(int i = 0; i<num_attribute; i++) {
				Attribute tmp_attribute = single_instance.attribute(i);
				if(tmp_attribute.isString()) {
					System.out.println(tmp_attribute.toString());
				}
			}

			
			single_instance.setClassIndex(single_instance.numAttributes() - 1);
			RandomForest rf = new RandomForest();
			
			weka.filters.supervised.attribute.AttributeSelection as = new  weka.filters.supervised.attribute.AttributeSelection();
		    Ranker ranker = new Ranker();
		    
		    
		    rf.buildClassifier(single_instance);
		    String info = rf.globalInfo();
		    System.out.println("Global info about Random Forest: "+info);
			
		    // cross validation step
		    int seed = 1;
		    int fold = 5;
		    Random rand = new Random(seed);   // create seeded number generator
		    Instances randData = new Instances(single_instance);   // create copy of original data
		    randData.randomize(rand);         // randomize data with number generator
		    
		    for (int n = 0; n < fold; n++) {
		    		Evaluation eval = new Evaluation(single_instance); 
		    		Instances train = randData.trainCV(fold, n, rand);
		    		
		    		Instances test = randData.testCV(fold, n);
		    		eval.evaluateModel(rf, test); 
		    	      // output fold statistics 
		    	    System.out.println("\nFold " + (n+1) + ":\n" + eval.toSummaryString()); 
		    		System.out.println(Integer.toString(n)+" --------------------------------------------------");
		    	
		    }
		    
		    
		    
		 // serialize model
		    ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(model_path));
		    oos.writeObject(rf);
		    oos.flush();
		    oos.close();
		    System.out.println(enzyme_name[s]+" model finished-------------------------------------------------");
		}
		
		
		
		
		
		
	}
	
	
	public static void main(String[] args) throws Exception {

		
		ArrayList<Instances> new_instance = create_phaseI_som_instance();
		create_som_model(new_instance,"load");
		
	}
}
