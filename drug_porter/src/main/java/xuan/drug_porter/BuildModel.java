package xuan.drug_porter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Random;

import org.apache.commons.lang3.ArrayUtils;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.fingerprint.SignatureFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;


import com.opencsv.CSVReader;
import com.opencsv.CSVWriter;

import weka.attributeSelection.AttributeSelection;
import weka.attributeSelection.BestFirst;
import weka.attributeSelection.CfsSubsetEval;
import weka.attributeSelection.Ranker;
import weka.classifiers.trees.RandomForest;
import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.ArffSaver;
import weka.core.converters.CSVLoader;
import weka.core.converters.ConverterUtils.DataSource;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.StringToNominal;
import xuan.drug_porter.descriptorUtils.*;

/**
 * csv file with all the data
 * 1. get csv file convert to instance (arff) format for weka
 * 1.1 should take the function that generate single smiles string to generate descriptor
 * 2. create couple of classifiers that take instance from step 1. give each function more control
 * 3. 
 * 
 * 
 * main: run multiple classifiers at same time can produce result without reduced feature
 * 		 run multiple classifiers at same time can produce result with reduced feature
 *       (optional) select the best algorithm to save the model for run classifier function
 *      
 * note: for each transporter, we need two models
 * 	     which biotransformer and cypreact miss => add inhibitor will compose the metabolism network
 * This class is only for build model purpose
 * 
 * 

 */
class BuildModel 
{
	
	
//	public static void test_fingerprint() throws CDKException {
//		SmilesParser temp_smiles = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//		IAtomContainer atom_container   = temp_smiles.parseSmiles("CCCCCCCC");
//		SignatureFingerprinter newfp = new SignatureFingerprinter();
//		BitSet newbit = newfp.getFingerprint(atom_container);
//		String bit = newbit.toString();
//	}
	//public Instances test_instances = null;
	
	public static String current_dir = System.getProperty("user.dir");
	private static String sep = File.pathSeparator;
	
	public static Instances get_training_instance(String csv_file, String model_type) throws Exception {
		
		
		
		String tempFile = current_dir+"/tmp_file/tmp_sdf.sdf";
		String instance_file_path = current_dir+"/tmp_file/tmp_instance_"+model_type+".arff";
		
		CSVReader reader = new CSVReader(new FileReader(csv_file));
	 	SDFWriter sdw  = new SDFWriter(new FileWriter(tempFile));
	 	
	 	
	 	String [] nextLine;
	     
	     //this loop will read all smile string, and convert it to sdf format of molecule
	     //then, write it back to sdf file SDFWriter sd
	 	
	     ArrayList<String> classification_class_list = new ArrayList<String>();
	     while ((nextLine = reader.readNext()) != null) {
	    	 String smile_string = nextLine[0];    //contain smile string
	    	 
	    	 if (!smile_string.isEmpty()) {
	    		 
	    		 classification_class_list.add(nextLine[1]);
	    		 SmilesParser temp_smiles = new SmilesParser(DefaultChemObjectBuilder.getInstance());
	    		 IAtomContainer atom_container   = temp_smiles.parseSmiles(smile_string);
	    		 StructureDiagramGenerator sdg = new StructureDiagramGenerator();
	    		 sdg.setMolecule(atom_container);
	    		 sdg.generateCoordinates();
	    		 IAtomContainer mole = sdg.getMolecule();
	    		 HashMap<Object,Object> properties = new HashMap<Object,Object>();
	    		 properties.put("SMILES", smile_string);
	    		 mole.addProperties(properties);
	 		
	    		 try {
	    			 sdw.write(mole);
	    		 } catch (Exception e) {
	    			 System.out.println(smile_string);
	    		 }
	    	 }  
	    	 else {
	    		 System.out.println(smile_string);
	    	 }
	 	}
	     sdw.close();
	     
	     
	     FeatureGeneration featureGeneration = new FeatureGeneration();
	     GenerateFeatureSingle GFS = new GenerateFeatureSingle();
	     GeneratingFeatures GF = new GeneratingFeatures();
	     
	     
	     
	     
	     IAtomContainerSet moleSet = featureGeneration.readFile(tempFile);
	     ArrayList<String[]> molecularFeatureList = new ArrayList<String[]>();
	     // generate feature name
	     ArrayList<Attribute> tmp_attributes = GF.generateAllAttributes(moleSet.getAtomContainer(0));
	     
	     FastVector<String> association = new FastVector<String>();
		 if (model_type.contains("substrate")) {
			System.out.println(model_type);
			association.addElement("substrate");
			association.addElement("non-substrate");
		 }
		 else if (model_type.contains("inhibitor")){
			 System.out.println(model_type);
			association.addElement("inhibitor");
			association.addElement("non-inhibitor");
		 }
		 Attribute class_attribute = new Attribute("Class",association);
		 tmp_attributes.add(class_attribute);
		 
		 Instances new_instance = new Instances("Rel",tmp_attributes,3000);
		 new_instance.setClassIndex(class_attribute.index()); 

	     
  
	     int num_of_attribute = tmp_attributes.size();
	     
	     FastVector feature_attribute = new FastVector(num_of_attribute);
	     for(Attribute a: tmp_attributes) {
	    	 	feature_attribute.addElement(a);
	    	 
	     }
	     
	     
	     
	     DecimalFormat df = new DecimalFormat("#.##");
	     for(int i = 0 ; i < moleSet.getAtomContainerCount(); i++) {
	    	 try {
	    		 
	    		 IAtomContainer mole = moleSet.getAtomContainer(i);	    		 
	    		 CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mole.getBuilder());
	    		 String Features = GFS.generateOneinstance(mole,"ALL");
	    		 
	    		 String class_  = classification_class_list.get(i);
	    		 Features = Features +","+class_;
	    		 
	    		 
		    	 String[] molecularFeature = Features.split(",");

	    		 Instance feature_set = new DenseInstance(num_of_attribute); // num_of_attribute=2281
	    		 
//	    		 System.out.println(molecularFeature.length);
	    		 for(int f=0; f < num_of_attribute; f++) {
	    			 Attribute tmp_attr = tmp_attributes.get(f);
	    			 if(tmp_attr.isNumeric()) {
	    				 feature_set.setValue(tmp_attr, Double.parseDouble(molecularFeature[f]));
	    			 }else if (tmp_attr.isNominal()) {
	    				 feature_set.setValue(tmp_attr, molecularFeature[f]);
	    			 }
	    			 
	    			 
	    		 }
	    		 
	    		 new_instance.add(feature_set);
	    	 }
	    	 catch (Exception e) {
	    		 System.out.println(e);
	    		 System.out.println(moleSet.getAtomContainer(i).getProperties());
	    	 }
	    		 
	     }
	     
	     System.out.println(new_instance.size());
	     File arff_instance = new File(instance_file_path);
		 ArffSaver arffSaver = new ArffSaver();
		 arffSaver.setInstances(new_instance);
		 arffSaver.setFile(arff_instance);
		 arffSaver.writeBatch();
			
		 reader.close();
		 return new_instance;
	     
		
		
	}
	
	
//	public static Instances get_test_instance(String csv_file) throws Exception {
//		String tempFile = current_dir+"\\tmp_file\\tmp_sdf.sdf";
//		String csv_writer_path = current_dir+"\\tmp_file\\tmp_csv.csv";
////		System.out.println(csv_writer_path);
////		System.out.println(tempFile);
//		
//		CSVReader reader = new CSVReader(new FileReader(csv_file));
//	 	SDFWriter sdw  = new SDFWriter(new FileWriter(tempFile));
//	 	CSVWriter writer = new CSVWriter(new FileWriter(csv_writer_path));
//	 	
//	 	
//	 	String [] nextLine;
//	     
//	     //this loop will read all smile string, and convert it to sdf format of molecule
//	     //then, write it back to sdf file SDFWriter sd
//	 	
//	     ArrayList<String> classification_class_list = new ArrayList<String>();
//	     while ((nextLine = reader.readNext()) != null) {
//	    	 String smile_string = nextLine[0];    //contain smile string
//	    	 classification_class_list.add(nextLine[1]);
// 	 		 SmilesParser temp_smiles = new SmilesParser(DefaultChemObjectBuilder.getInstance());
// 	 		 IAtomContainer atom_container   = temp_smiles.parseSmiles(smile_string);
//	 		 StructureDiagramGenerator sdg = new StructureDiagramGenerator();
//	 		 sdg.setMolecule(atom_container);
//	 		 sdg.generateCoordinates();
//	 		 IAtomContainer mole = sdg.getMolecule();
//	 		 HashMap<Object,Object> properties = new HashMap<Object,Object>();
//	 		 properties.put("SMILES", smile_string);
//	 		 mole.addProperties(properties);
//	 		
//	 		 try {
//	 			 sdw.write(mole);
//	 		 } catch (Exception e) {
//	 			 System.out.println(smile_string);
//	 		 }
//	 	}
//	     sdw.close();
//	     
//	     
//	     String Attributes = null;
//	     
//	     FeatureGeneration featureGeneration = new FeatureGeneration();
//	     GenerateFeatureSingle GFS = new GenerateFeatureSingle();
//	     GeneratingFeatures GF = new GeneratingFeatures();
//	     
//	     
//	     
//	     
//	     IAtomContainerSet moleSet = featureGeneration.readFile(tempFile);
//	     Attributes = GF.generateAllFeatures(moleSet.getAtomContainer(0),"molecularFeatures");
//	     Attributes = Attributes +"Class"; // add the classification class
//	     String[] attribute_list = Attributes.split(",");
//	     writer.writeNext(attribute_list);
//
//	     //String attribute = GFS.generate_all_cdk_molecular_descriptor(moleSet.getAtomContainer(0), feature_type,"name");
//	     // Iterate through all the molecule from AtomContainer
//	     for(int i = 0 ; i < moleSet.getAtomContainerCount(); i++) {
//	    	
//			IAtomContainer mole = moleSet.getAtomContainer(i);
//			System.out.println(mole.getProperties());
//			
//			
//			CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mole.getBuilder());
//	    	
//			String Features = GFS.generateOneinstance(mole,"ALL");
//			
//			// GFS.generate_all_cdk_molecular_descriptor will return all molecular feature from cdk library
//			//String Features_all = GFS.generate_all_cdk_molecular_descriptor(mole, feature_type,"value")
//			Features = Features +",?";
//			String[] molecularFeature = Features.split(",");
//
//			writer.writeNext(molecularFeature);
//		}
//	    
//	    writer.close();
//	    
//	    
//		File checkFile = new File(tempFile);
//		if(checkFile.exists()) {
//			checkFile.delete();
//			System.out.println("Temp File deleted");
//		}
//		reader.close();
//		// read csv -> convert to arff
//		// load CSV
//		CSVLoader loader = new CSVLoader();
//		loader.setSource(new File(csv_writer_path));
//		Instances data = loader.getDataSet();	
//		
//		File checkFile2 = new File(csv_writer_path);
//		if(checkFile2.exists()) {
//			checkFile2.delete();
//			System.out.println("temp_csv.csv File deleted");    
//		}
//		
//		
//
//		return data;
//	}
	
	
	
	/*
	 * classifier function
	 */
	
	// Random Forest
	public static void run_random_forest(Instances training_instance, Instances test_instances,String model_type) throws Exception {
		
		
		String model_path =  current_dir+"/tmp_file/RandomForest_"+model_type+".model";
        
		int num_attribute = training_instance.numAttributes();
		for(int i = 0; i<num_attribute; i++) {
			Attribute tmp_attribute = training_instance.attribute(i);
			if(tmp_attribute.isString()) {
				System.out.println(tmp_attribute.toString());
			}
		}

		
		training_instance.setClassIndex(training_instance.numAttributes() - 1);
		RandomForest rf = new RandomForest();
		
		weka.filters.supervised.attribute.AttributeSelection as = new  weka.filters.supervised.attribute.AttributeSelection();
	    Ranker ranker = new Ranker();
	    
	    
	    rf.buildClassifier(training_instance);
	    String info = rf.globalInfo();
	    System.out.println("Global info about Random Forest: "+info);
		
	    // cross validation step
	    int seed = 1;
	    int fold = 5;
	    Random rand = new Random(seed);   // create seeded number generator
	    Instances randData = new Instances(training_instance);   // create copy of original data
	    randData.randomize(rand);         // randomize data with number generator
	    
	    for (int n = 0; n < fold; n++) {
	    		Instances train = randData.trainCV(fold, n, rand);
	    		
	    		Instances test = randData.testCV(fold, n);
	    		System.out.println(Integer.toString(n)+" run-------------------------------------------------");
	    	
	    }
	    
	    
	    
	 // serialize model
	    ObjectOutputStream oos = new ObjectOutputStream(
	                               new FileOutputStream(model_path));
	    oos.writeObject(rf);
	    oos.flush();
	    oos.close();
		
	}
	
	
	// Support Vector Machine
	public void run_SVM(Instances training_instance, Instances test_instances) {
		

		
	}
	
	
	// Logistic Regression
	public void run_logistic_regression(Instances training_instance, Instances test_instances) {
		

		
	}
	
	
	// NN
	public void run_NN(Instances traning_instance,Instances test_instances) {
		
	}
	
	
	// navie_bayes
	public void run_navie_bayes(Instances traning_instance,Instances test_instances) {
		
	}
	
	
	
	// reduce feature
	static Instances performFeatureExtraction(Instances data) {
	    System.out.println(data.numAttributes());
	    CfsSubsetEval evaluator = new CfsSubsetEval();
	    int k = data.numAttributes();

	    BestFirst ranker = new BestFirst(); //new Ranker();

	    try {
	      AttributeSelection selector = new AttributeSelection();
	      selector.setSearch(ranker);
	      selector.setEvaluator(evaluator);
	      selector.SelectAttributes(data);

	      // Transform data into eigenvector basis.
	      Instances transformedData = selector.reduceDimensionality(data);
	      System.out.println(transformedData.numAttributes());
	      return transformedData;
	    } catch (Exception e) {
	      e.printStackTrace();
	      return data;
	    }

	}
	
	
    public static void main( String[] args ) throws Exception
    {
//    	String csv_path = "C:\\Users\\Danis\\Desktop\\inhibitor.csv";
//    	String csv_path_test = "C:\\Users\\Danis\\Desktop\\test_tmp_csv.csv";
    	 
    		String csv_path   = "/Users/xuan/Desktop/tmp_transporter/substrate.csv";
		String[] splits = csv_path.split("/");
		String model_type = splits[splits.length-1];
		model_type = model_type.replace(".csv","");
		Instances training_instance = get_training_instance(csv_path,model_type);
		//    Instances test_instance = get_test_instance(csv_path_test);
		run_random_forest(training_instance,training_instance,model_type);
		
		// inhibitor model
		String csv_path_inhibitor   = "/Users/xuan/Desktop/tmp_transporter/inhibitor.csv";
		String[] splits_inhibitor = csv_path_inhibitor.split("/");
		String model_type_2 = splits_inhibitor[splits_inhibitor.length-1];
		model_type_2 = model_type_2.replace(".csv","");
		Instances training_instance_inhibitor = get_training_instance(csv_path_inhibitor,model_type_2);
		//    Instances test_instance = get_test_instance(csv_path_test);
		run_random_forest(training_instance_inhibitor,training_instance_inhibitor,model_type_2);
		
		
		
		
		
//    		String instance_file_path = current_dir+"/tmp_file/tmp_instance.arff";
//    		File arff_instance = new File(instance_file_path);
//    	if (!arff_instance.exists())
//        {
//    		
//    			String csv_path   = args[0];
//    			String model_type = args[1];
//    		
//    			Instances training_instance = get_training_instance(csv_path,model_type);
////            Instances test_instance = get_test_instance(csv_path_test);
//            run_random_forest(training_instance,training_instance);
//        }
//    	else {
//    		DataSource source = new DataSource(instance_file_path);
//    		Instances data = source.getDataSet();
//    		run_random_forest(data,data);
//    	}
        
        
    }
}
