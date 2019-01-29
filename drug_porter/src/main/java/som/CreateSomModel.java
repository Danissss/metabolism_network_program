package som;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;

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

	

	public static String current_dir = System.getProperty("user.dir");
	
	public static Instances create_phaseI_som_instance() throws CDKException, IOException {
		
//		ArrayList<ArrayList<ArrayList<String>>> dataset = create_phaseI_som_3d_arrayList();
		ArrayList<ArrayList<HashMap<IAtomContainer,ArrayList<ArrayList<String>>>>> dataset =  create_phaseI_som_arraylist();
		int nearest_atom = 3;	// find nearest 4 atoms
		// previous step could be done by just call create_phaseI_som3d_arraylist() function
		// create instance;
		String[] enzyme_name = new String[] {"1A2","2A6","2B6","2C8","2C9","2C19","2D6","2E1","3A4"}; 
		
		
		ArrayList<Attribute> attribute_name = generate_attribute_name(nearest_atom+1, 29); // return 4*29=116 attribute name;
	    
		 
		for(int i = 0; i< dataset.size(); i++) {
			String arff_path = current_dir+"/tmp_file/"+enzyme_name[i]+"_som.arff";
			
			
			// create the instance and attribute
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
					int num_atom = mole.getAtomCount();
					ArrayList<ArrayList<String>> all_nearest_atom_set = GetAtomicDescriptors.getNearestAtoms(mole);  // return nearest atom list in atom order
					ArrayList<ArrayList<String>> instances = entry.getValue();
					for(int k = 0; k<instances.size(); k++) {
						
						// iterate each instance;
						// convert each instance to weka instance;
						int som = Integer.parseInt(instances.get(k).get(0));
						String is_som = instances.get(k).get(1);
						ArrayList<String> nearest_atom_set = all_nearest_atom_set.get(som-1);   // get index of atom but need to minus 1 because of chemsketch
						List<IAtom> atoms_set = new ArrayList<IAtom>();
						
						// add the # nearest atom into atom list for extract the atom descriptor;
						atoms_set.add(mole.getAtom(som));  									// add the original atoms
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
			    		 			feature_set.setValue(tmp_attr, single_instance_value.get(f));
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
			
			
			System.out.println(arff_path+" finished!======================!");
		}
		
		return null;
		
	}
	/**
	 * 
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
		
		
		ArrayList<ArrayList<String>> list_1A2 = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> list_2A6 = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> list_2B6 = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> list_2C8 = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> list_2C9 = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> list_2C19 = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> list_2D6 = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> list_2E1 = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> list_3A4 = new ArrayList<ArrayList<String>>();
		
		
		// hashmap<mol,ArrayList<ArrayList<>>
		
		
		
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
					array_1A2.add(Integer.toString(value_int));
					array_1A2.add("Yes");
					list_1A2.add(array_1A2);		// add to the list_1A2 (ArrayList<ArrayList<String>>)
					array_1A2.clear(); 			// clean the array_1A2
					flag_1A2 = 1;
				}
				else if (key.contains("2A6")) {
					int value_int = Integer.valueOf(entry.getValue().toString());
					atom_list_2A6.remove(value_int-1);
					array_2A6.add(Integer.toString(value_int));
					array_2A6.add("Yes");
					list_2A6.add(array_2A6);
					array_2A6.clear();
					flag_2A6 = 1;
				}
				else if (key.contains("2B6")) {
					int value_int = Integer.valueOf(entry.getValue().toString());
					atom_list_2B6.remove(value_int-1);
					array_2B6.add(Integer.toString(value_int));
					array_2B6.add("Yes");
					list_2B6.add(array_2B6);
					array_2B6.clear();
					flag_2B6 = 1;
				}
				else if (key.contains("2C8")) {
					int value_int = Integer.valueOf(entry.getValue().toString());
					atom_list_2C8.remove(value_int-1);
					array_2C8.add(Integer.toString(value_int));
					array_2C8.add("Yes");
					list_2C8.add(array_2C8);
					array_2C8.clear();
					flag_2C8 = 1;
				}
				else if (key.contains("2C9")) {
					int value_int = Integer.valueOf(entry.getValue().toString());
					atom_list_2C9.remove(value_int-1);
					array_2C9.add(Integer.toString(value_int));
					array_2C9.add("Yes");
					list_2C9.add(array_2C9);
					array_2C9.clear();
					flag_2C9 = 1;
				}
				else if (key.contains("2C19")) {
					int value_int = Integer.valueOf(entry.getValue().toString());
					atom_list_2C19.remove(value_int-1);
					array_2C19.add(Integer.toString(value_int));
					array_2C19.add("Yes");
					list_2C19.add(array_2C19);
					array_2C19.clear();
					flag_2C19 = 1;
				}
				else if (key.contains("2D6")) {
					int value_int = Integer.valueOf(entry.getValue().toString());
					atom_list_2D6.remove(value_int-1);
					array_2D6.add(Integer.toString(value_int));
					array_2D6.add("Yes");
					list_2D6.add(array_2D6);
					array_2D6.clear();
					flag_2D6 = 1;
				}
				else if (key.contains("2E1")) {
					int value_int = Integer.valueOf(entry.getValue().toString());
					atom_list_2E1.remove(value_int-1);
					array_2E1.add(Integer.toString(value_int));
					array_2E1.add("Yes");
					list_2E1.add(array_2E1);
					array_2E1.clear();
					flag_2E1 = 1;
				}
				else if (key.contains("3A4")) {
					int value_int = Integer.valueOf(entry.getValue().toString());
					atom_list_3A4.remove(value_int-1);
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
					array_1A2.add(Integer.toString(atom_list_1A2.get(left_atom)));
					array_1A2.add("No");
					list_1A2.add(array_1A2);
					array_1A2.clear();
					
				}
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
				properties.put(mol,list_1A2);
				hash_1A2.add(properties);
				
				
			}
			
			if(flag_2A6 == 1) {
				ArrayList<String> array_2A6 = new ArrayList<String>();
				for (int left_atom = 0; left_atom < atom_list_2A6.size();left_atom++) {
					array_2A6.add(Integer.toString(atom_list_2A6.get(left_atom)));
					array_2A6.add("No");
					list_2A6.add(array_2A6);
					array_2A6.clear();
				}
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
				properties.put(mol,list_2A6);
				hash_2A6.add(properties);
			}
			
			if(flag_2B6 == 1) {
				ArrayList<String> array_2B6 = new ArrayList<String>();
				for (int left_atom = 0; left_atom < atom_list_2B6.size();left_atom++) {
					array_2B6.add(Integer.toString(atom_list_2B6.get(left_atom)));
					array_2B6.add("No");
					list_2B6.add(array_2B6);
				}
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
				properties.put(mol,list_2B6);
				hash_2B6.add(properties);
			}
			
			if(flag_2C8 == 1) {
				ArrayList<String> array_2C8 = new ArrayList<String>();
				for (int left_atom = 0; left_atom < atom_list_2C8.size();left_atom++) {	
					array_2C8.add(Integer.toString(atom_list_2C8.get(left_atom)));
					array_2C8.add("No");
					list_2C8.add(array_2C8);
					array_2C8.clear();
				}
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
				properties.put(mol,list_2C8);
				hash_2C8.add(properties);
			}
			
			if(flag_2C9 == 1) {
				ArrayList<String> array_2C9 = new ArrayList<String>();
				for (int left_atom = 0; left_atom < atom_list_2C9.size();left_atom++) {
					array_2C9.add(Integer.toString(atom_list_2C9.get(left_atom)));
					array_2C9.add("No");
					list_2C9.add(array_2C9);
					array_2C9.clear();
				}
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
				properties.put(mol,list_2C9);
				hash_2C9.add(properties);
			}
			
			if(flag_2C19 == 1) {
				ArrayList<String> array_2C19 = new ArrayList<String>();
				for (int left_atom = 0; left_atom < atom_list_2C19.size();left_atom++) {
					array_2C19.add(Integer.toString(atom_list_2C19.get(left_atom)));
					array_2C19.add("No");
					list_2C19.add(array_2C19);
					array_2C19.clear();
				}
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
				properties.put(mol,list_2C19);
				hash_2C19.add(properties);
			}
			
			if(flag_2D6 == 1) {
				ArrayList<String> array_2D6 = new ArrayList<String>();
				for (int left_atom = 0; left_atom < atom_list_2D6.size();left_atom++) {
					array_2D6.add(Integer.toString(atom_list_2D6.get(left_atom)));
					array_2D6.add("No");
					list_2D6.add(array_2D6);
					array_2D6.clear();
				}
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
				properties.put(mol,list_2D6);
				hash_2D6.add(properties);
			}
			
			if(flag_2E1 == 1) {
				ArrayList<String> array_2E1 = new ArrayList<String>();
				for (int left_atom = 0; left_atom < atom_list_2E1.size();left_atom++) {
					array_2E1.add(Integer.toString(atom_list_2E1.get(left_atom)));
					array_2E1.add("No");
					list_2E1.add(array_2E1);
					array_2E1.clear();
				}
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
				properties.put(mol,list_2E1);
				hash_2E1.add(properties);
			}
			
			if(flag_3A4 == 1) {
				ArrayList<String> array_3A4 = new ArrayList<String>();
				for (int left_atom = 0; left_atom < atom_list_3A4.size();left_atom++) {
					array_3A4.add(Integer.toString(atom_list_3A4.get(left_atom)));
					array_3A4.add("No");
					list_3A4.add(array_3A4);
					array_3A4.clear();
				}
				
				HashMap<IAtomContainer,ArrayList<ArrayList<String>>> properties = new HashMap<IAtomContainer,ArrayList<ArrayList<String>>>();
				properties.put(mol,list_3A4);
				hash_3A4.add(properties);
			}
			
			list_1A2.clear();
			list_2A6.clear();
			list_2B6.clear();
			list_2C8.clear();
			list_2C9.clear();
			list_2C19.clear();
			list_2D6.clear();
			list_2E1.clear();
			list_3A4.clear();
		}
		
		
		
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
	
	public static void main(String[] args) throws FileNotFoundException, CDKException {
		
//		create_phaseI_som_instance();
		
	}
}
