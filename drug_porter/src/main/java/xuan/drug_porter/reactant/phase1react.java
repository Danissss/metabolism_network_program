package xuan.drug_porter.reactant;

/** 
 * Authors: Siyang Tian
 * Class Description:
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;


import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import weka.classifiers.Classifier;
import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;


import xuan.drug_porter.descriptorUtils.SdfToSample;

public class phase1react {
	
	
	protected static String current_dir = System.getProperty("user.dir");
	private String sep = File.pathSeparator;
	
	public static void main(String[] args) throws Exception{
		if(args.length<4){
			System.out.println("Don't have enough arguments");
			return;
		}
		else if(args.length>4){
			System.out.println("# of arguments is more than 4");
			return;			
		}
		String supportFoldPath = args[0]; //Path of the folder of the supportfiles ( .model and supportfiles.csv)
		String input = args[1]; //The path of the input .sdf/.csv file
		String output = args[2]; //The path of the ouput file, it can be either a .sdf or a .csv file
		String cyp = args[3]; //The target CYPs

		
		
		
	}
	public void makePrediction(String supportFoldPath, String input, String output, String cyp) throws Exception {
		
		
		String[] cypList = cyp.split(","); 
		SdfToSample sf = new SdfToSample();
		phase1react test = new phase1react();
		IAtomContainerSet inputMolecules = sf.createIAtomContainerSet(input);
		SdfToSample sdfTool = new SdfToSample();		
		Instances testSet = sdfTool.generateAllInstances(inputMolecules);
		String supportfile = "Determined by the input";
		//The predictedResult ArrayList is used to store all 9 predicted results;
		ArrayList<HashMap<String,String>> predictedResult = new ArrayList<HashMap<String,String>>();
		
		//Create arrayList predictedResult with size = inputMolecules.getAtomContainerCount(), store results fore every molecule
		predictedResult = test.initPreResults(predictedResult,inputMolecules.getAtomContainerCount());
		
		String model = "To be choosen";
		
		// predictResult is ArrayList<HashMap<String,String>>
		if(cypList.length>1){
			
			predictedResult = test.makeMultiPrediction(supportFoldPath, cyp, testSet, predictedResult);
		}
		
		else if(cyp.contains("1A2")){	
			model = supportFoldPath + "supportfiles/CYP1A2/model/1A2_NR.model";
			supportfile = supportFoldPath + "supportfiles/CYP1A2/supportfile.csv";
			model = current_dir+"/Model/Phase1/1A2_NR.model";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makePrediction("1A2",model,testSet,supportfile,predictedResult );
		}
		else if(cyp.contains("3A4")){
			model = supportFoldPath + "supportfiles/CYP3A4/model/3A4_NR.model";
		
			supportfile = supportFoldPath + "supportfiles/CYP3A4/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makePrediction("3A4",model,testSet,supportfile,predictedResult );
		}
		else if(cyp.contains("2C9")){			
			supportfile = supportFoldPath + "supportfiles/CYP2C9/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makeEnsemblePrediction("2C9",testSet,supportFoldPath,supportfile,predictedResult );
		}
		else if(cyp.contains("2C19")){
			model = supportFoldPath + "supportfiles/CYP2C19/model/2C19_NR.model";
			supportfile = supportFoldPath + "supportfiles/CYP2C19/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makePrediction("2C19",model,testSet,supportfile,predictedResult );
		}
		else if(cyp.contains("2E1")){
			model = supportFoldPath + "supportfiles/CYP2E1/model/2E1_NR.model";			
			supportfile = supportFoldPath + "supportfiles/CYP2E1/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makePrediction("2E1",model,testSet,supportfile,predictedResult );
		}
		else if(cyp.contains("2D6")){
			supportfile = supportFoldPath + "supportfiles/CYP2D6/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makeEnsemblePrediction("2D6",testSet,supportFoldPath,supportfile,predictedResult );
		}
		else if(cyp.contains("2A6")){
			model = supportFoldPath + "supportfiles/CYP2A6/model/2A6_NR.model";
			supportfile = supportFoldPath + "supportfiles/CYP2A6/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makePrediction("2A6",model,testSet,supportfile,predictedResult );
		}
		else if(cyp.contains("2B6")){
			model = supportFoldPath + "supportfiles/CYP2B6/model/2B6_NR.model";
			supportfile = supportFoldPath + "supportfiles/CYP2B6/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makePrediction("2B6",model,testSet,supportfile,predictedResult);
		}
		else if(cyp.contains("2C8")){
			model = supportFoldPath + "supportfiles/CYP2C8/model/2C8_NR.model";			
			supportfile = supportFoldPath + "supportfiles/CYP2C8/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makePrediction("2C8",model,testSet,supportfile,predictedResult );
		}
		else{
			System.out.println("Be added soon");
			return;
		}
		if(output.contains("sdf")){
			test.outputResultIAtomContainerSet(inputMolecules, output, predictedResult);
		}
		else if(output.contains("csv")){
			test.outPutCsv(inputMolecules, output, predictedResult);
		}
		else System.out.println("Cannot detect the output file format");
	}
	
	/**
	 * Predict whether the given test molecules in the sdf/csv file are reactants. This is called when more than one Cyps are input
	 * @param supportFoldPath
	 * @param cyp the input CYPs, spited by "," 
	 * @param inputMolecules: The input compounds
	 * @param predictedResult: The predicted results for all input compounds.
	 * @return updated predictedResult
	 * @throws Exception
	 */
	public ArrayList<HashMap<String,String>> makeMultiPrediction(String supportFoldPath, String cyp, Instances testSet, ArrayList<HashMap<String,String>> predictedResult) throws Exception{
		ArrayList<HashMap<String,String>> resultList = predictedResult;
		if(cyp.contains("1A2")){
			String model = supportFoldPath + "supportfiles/CYP1A2/model/1A2_NR.model";
			String supportfile = supportFoldPath + "supportfiles/CYP1A2/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			resultList = makePrediction("1A2",model, testSet,supportfile,resultList );
		}
		if(cyp.contains("3A4")){
			String model = supportFoldPath + "supportfiles/CYP3A4/model/3A4_NR.model";
			String supportfile = supportFoldPath + "supportfiles/CYP3A4/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			resultList = makePrediction("3A4",model,testSet,supportfile,resultList );
		}
		if(cyp.contains("2B6")){
			String model = supportFoldPath + "supportfiles/CYP2B6/model/2B6_NR.model";
			String supportfile = supportFoldPath + "supportfiles/CYP2B6/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			resultList = makePrediction("2B6",model,testSet,supportfile,resultList );
		}
		if(cyp.contains("2E1")){
			String model = supportFoldPath + "supportfiles/CYP2E1/model/2E1_NR.model";
			String supportfile = supportFoldPath + "supportfiles/CYP2E1/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			resultList = makePrediction("2E1",model,testSet,supportfile,resultList );
		}
		if(cyp.contains("2C9")){
			String supportfile = supportFoldPath + "supportfiles/CYP2C9/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			resultList = makeEnsemblePrediction("2C9",testSet,supportFoldPath,supportfile,predictedResult );
		}
		if(cyp.contains("2C19")){
			String model = supportFoldPath + "supportfiles/CYP2C19/model/2C19_NR.model";
			String supportfile = supportFoldPath + "supportfiles/CYP2C19/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			resultList = makePrediction("2C19",model,testSet,supportfile,resultList );
		}
		if(cyp.contains("2D6")){
			String supportfile = supportFoldPath + "supportfiles/CYP2D6/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			resultList = makeEnsemblePrediction("2D6",testSet, supportFoldPath,supportfile,predictedResult);
		}
		if(cyp.contains("2C8")){
			String model = supportFoldPath + "supportfiles/CYP2C8/model/2C8_NR.model";
			String supportfile = supportFoldPath + "supportfiles/CYP2C8/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			resultList = makePrediction("2C8",model,testSet,supportfile,resultList );
		}
		if(cyp.contains("2A6")){
			String model = supportFoldPath + "supportfiles/CYP2A6/model/2A6_NR.model";
			String supportfile = supportFoldPath +  "supportfiles/CYP2A6/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			resultList = makePrediction("2A6",model,testSet,supportfile,resultList );
		}
		return resultList;
	}
	
	/**
	 * Predict whether the given molecules in the sdf/csv file are reactant or not for one CYP by using ensemble method
	 * @param cyp
	 * @param inputMolecules
	 * @param supportfile: The path to the supportfiles, including .models and supportfile.csv
	 * @param predictedResult
	 * @return updated predictedResult
	 * @throws Exception
	 */
	public ArrayList<HashMap<String,String>> makeEnsemblePrediction(String cyp, Instances testSet, String supportFoldPath, String supportfile, ArrayList<HashMap<String,String>> predictedResult) throws Exception{
		
		ArrayList<HashMap<String,String>> resultList = predictedResult;
		ArrayList<Classifier> enSembles = new ArrayList<Classifier>();
		//Load the ensemble models
		if(cyp == "2D6"){
			for(int i = 0; i < 5; i++){
				String model = supportFoldPath + "supportfiles/CYP" + cyp + "/model/" + cyp + "_NR_" + i +".model";
				Classifier cls = (Classifier) weka.core.SerializationHelper.read(model);
				enSembles.add(cls);
			}
		
		}
		if(cyp == "2C9"){
			for(int i = 0; i < 6; i++){
				String model = supportFoldPath + "supportfiles/CYP" + cyp + "/model/" + cyp + "_NR_" + i +".model";
				Classifier cls = (Classifier) weka.core.SerializationHelper.read(model);
				enSembles.add(cls);
			}
		
		}
		
		File supfile = new File(supportfile);
		FileReader spfr = new FileReader(supfile);
		BufferedReader spbr = new BufferedReader(spfr);
		
		//Declare attributes arraylist,mean arraylist, max arraylist, min arraylist, will be used in normalization
		ArrayList<String> attList = new ArrayList<String>();
		ArrayList<String> meanList = new ArrayList<String>();
		ArrayList<String> maxList = new ArrayList<String>();
		ArrayList<String> minList = new ArrayList<String>();
		int counter = 0;
		//update mean, max, min values of each Attribute.
		String supLine = spbr.readLine();//Skip the first line that contains all the titles
		while((supLine = spbr.readLine())!=null){
			String[] elements = supLine.split(",");
			attList.add(elements[0]);
			meanList.add(elements[1]);
			maxList.add(elements[2]);
			minList.add(elements[3]);
		}
		int countR = 0;
		int countN = 0;
		Instances matched = matchAttributes(testSet,attList,meanList,maxList,minList);
		matched.setClassIndex(matched.numAttributes()-1);
		System.out.println("--------------------------------------");
		for(int i = 0; i<matched.size();i++){
			counter++;
			Instance oneSample = matched.get(i);
			double tempResult = 0.0;
			for(int k = 0; k < enSembles.size(); k++){
				Classifier cls = enSembles.get(k);
				tempResult = tempResult + cls.classifyInstance(oneSample);
			}
			double result = 0.0;
			if(tempResult/enSembles.size() >= 0.5){
				result = 1.0;
			}
			else result = 0.0;
			String pred = "";
		
			if(result==0.0){
				countN++;
				System.out.println(cyp +", " +"Mole" + (i+1) + ": N");
				pred = "N";
			}
			else if(result == 1.0){
				countR++;
				System.out.println(cyp +", " +"Mole" + (i+1) + ": R");
				pred = "R";
			}
							
			resultList.get(i).replace(cyp,pred);//Update the predicted result of the corresponding CYP
		}
		//System.out.println(cyp + " " + "N:" + countN + " R:" + countR );
		return resultList;
	}	
	
	
	/**
	 * Predict whether the given molecules in the sdf/csv file are reactants or not for one CYP by using non-ensemble methods
	 * @param cyp
	 * @param model
	 * @param testSet
	 * @param supportfile
	 * @param predictedResult
	 * @return updated PredictedResults
	 * @throws Exception
	 */

	public ArrayList<HashMap<String,String>> makePrediction(String cyp, String model, Instances testSet, String supportfile, ArrayList<HashMap<String,String>> predictedResult) throws Exception{
		
		ArrayList<HashMap<String,String>> resultList = predictedResult;			
		Classifier cls = (Classifier) weka.core.SerializationHelper.read(model);//Load the model
				
		File supfile = new File(supportfile);
		FileReader spfr = new FileReader(supfile);
		BufferedReader spbr = new BufferedReader(spfr);
		
		//Create attribute arraylist,mean arraylist, max arraylist, min arraylist
		ArrayList<String> attList = new ArrayList<String>();
		ArrayList<String> meanList = new ArrayList<String>();
		ArrayList<String> maxList = new ArrayList<String>();
		ArrayList<String> minList = new ArrayList<String>();
		//Write the output 
		int counter = 0;

		//attribute	mean	max	min. Skip the first line that contains all the titles
		String supLine = spbr.readLine();
		while((supLine = spbr.readLine())!=null){
			String[] elements = supLine.split(",");
			attList.add(elements[0]);
			meanList.add(elements[1]);
			maxList.add(elements[2]);
			minList.add(elements[3]);
		}
		int countR = 0;
		int countT = 0;
		Instances matched = matchAttributes(testSet,attList,meanList,maxList,minList);
		matched.setClassIndex(matched.numAttributes()-1);
		System.out.println("--------------------------------------");
		for(int i = 0; i<matched.size();i++){
			counter++;
			Instance oneSample = matched.get(i);
			double result = cls.classifyInstance(oneSample);
			String pred = "";
		
			if(result==0.0){
				countT++;
				System.out.println(cyp +", " +"Mole" + (i+1) +": N");
				pred = "N";
			}
			else if(result == 1.0){
				countR++;
				System.out.println(cyp +", " +"Mole" + (i+1) + ": R");
				pred = "R";
			}			
			
			resultList.get(i).replace(cyp,pred);
		}
		//System.out.println("N:" + countT + " R:" + countR );
		return resultList;
	}
	
	/**
	 * Write input molecules with predicted results into the target sdf file
	 * @param rawMolecules: The original input compounds
	 * @param outSdfPath: The path to the target sdf file
	 * @param predictedResult: predicted results for all compounds
	 * @throws Exception
	 */

	public void outputResultIAtomContainerSet(IAtomContainerSet rawMolecules, String outSdfPath, ArrayList<HashMap<String,String>> predictedResult) throws Exception{
		IAtomContainerSet resultMole = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		String sdfOutput = outSdfPath;
		SDFWriter sdfWriter  = new SDFWriter(new FileWriter(sdfOutput)); 
		
		IAtomContainerSet resultIAtomContainerSet = getResultIAtomContainerSet(rawMolecules, predictedResult);
		for(int i = 0; i < resultIAtomContainerSet.getAtomContainerCount(); i++){			
			IAtomContainer outMole = resultIAtomContainerSet.getAtomContainer(i);
			sdfWriter.write(AtomContainerManipulator.removeHydrogens(outMole));
		}
		sdfWriter.close();
	}
	
	/**
	 * Return IAtomContainerSet that contains all predicted results and molecules properties for the input molecules
	 * @param rawMolecules: The original input compounds
	 * @param predictedResult: predicted results for all compounds
	 * @return input compounds with predicted results
	 * @throws Exception
	 */

	public IAtomContainerSet getResultIAtomContainerSet(IAtomContainerSet rawMolecules, ArrayList<HashMap<String,String>> predictedResult) throws Exception{
		IAtomContainerSet resultMole = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		for(int i = 0; i < rawMolecules.getAtomContainerCount(); i++){
			IAtomContainer oneMole = rawMolecules.getAtomContainer(i);
			IAtomContainer outMole = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
			outMole = oneMole.clone();
			outMole.setProperties(null);
			//Create proper list
			Map<Object, Object> CypPre = new LinkedHashMap<Object, Object>();
			//Add inchikey
			// Generate factory - throws CDKException if native code does not load
			InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
			// Get InChIGenerator
			InChIGenerator gen1 = factory.getInChIGenerator(oneMole);
			String inchikey1 = gen1.getInchiKey();
			CypPre.put("InchiKey", inchikey1);
			
			//Add SMILEs
			SmilesGenerator sg = new SmilesGenerator().unique();
			String smile1 = sg.create(oneMole);
			CypPre.put("SMILES", smile1);
			
			//Add Title
			String title1 = oneMole.getProperty(CDKConstants.TITLE); 
			if(title1 != null && !title1.isEmpty()){
				CypPre.put("cdk:Title", title1);
			}
			else {
				title1 = "molecule:" + (i+1);
				CypPre.put("cdk:Title", title1);
			}			
			CypPre.put("1A2", predictedResult.get(i).get("1A2"));
			CypPre.put("2B6", predictedResult.get(i).get("2B6"));
			CypPre.put("2A6", predictedResult.get(i).get("2A6"));
			CypPre.put("2C8", predictedResult.get(i).get("2C8"));
			CypPre.put("2C9", predictedResult.get(i).get("2C9"));
			CypPre.put("2C19",predictedResult.get(i).get("2C19"));
			CypPre.put("2D6", predictedResult.get(i).get("2D6"));
			CypPre.put("2E1", predictedResult.get(i).get("2E1"));
			CypPre.put("3A4", predictedResult.get(i).get("3A4"));

			Map<Object, Object> oldProperties = oneMole.getProperties();

			//Filter out the null mappings from the oldProperties
			HashMap<Object,Object> filteredProperties = new HashMap<Object,Object>();
			for(Object key : oldProperties.keySet()){
				if(key.equals("2C9")&&(oldProperties.get(key)!=null)){
					CypPre.put("2C9", oldProperties.get(key));		
				}
				if(key.equals("1A2")&&(oldProperties.get(key)!=null)){
					CypPre.put("1A2", oldProperties.get(key));		
				}
				if(key.equals("2A6")&&(oldProperties.get(key)!=null)){
					CypPre.put("2A6", oldProperties.get(key));		
				}
				if(key.equals("2B6")&&(oldProperties.get(key)!=null)){
					CypPre.put("2B6", oldProperties.get(key));		
				}
				if(key.equals("2C8")&&(oldProperties.get(key)!=null)){
					CypPre.put("2C8", oldProperties.get(key));		
				}
				if(key.equals("2C19")&&(oldProperties.get(key)!=null)){
					CypPre.put("2C19", oldProperties.get(key));		
				}
				if(key.equals("2E1")&&(oldProperties.get(key)!=null)){
					CypPre.put("2E1", oldProperties.get(key));		
				}
				if(key.equals("3A4")&&(oldProperties.get(key)!=null)){
					CypPre.put("3A4", oldProperties.get(key));		
				}
				if(key.equals("2D6")&&(oldProperties.get(key)!=null)){
					CypPre.put("2D6", oldProperties.get(key));		
				}
				else if(oldProperties.get(key)!=null&&(!oldProperties.get(key).equals("null"))){
					filteredProperties.put(key, oldProperties.get(key));
				}
			}
			oldProperties = filteredProperties;
			Set<Object> keySets = oldProperties.keySet();
			// The section below checks whether title, Smiles and InchiKey from the origin sdf file are the same as the one generated by our tool
			if(oldProperties.keySet().contains("cdk:Title")){
				String checkTitle = (String) oldProperties.get("cdk:Title");
				if(!checkTitle.equals(title1)){
					System.out.println("The " + i + "th mole's title is different from the title generated by our tool");
					System.out.println("Origin title: " + checkTitle);
					System.out.println("OurGen title: " + title1);
					//oldProperties.remove("cdk:Title", null);
					oldProperties.put("cdk:Title",title1);
					
				}
			}
			if(oldProperties.keySet().contains("InchiKey")){
				String checkInchikey = (String) oldProperties.get("InchiKey");
				if(!checkInchikey.equals(inchikey1)){
					System.out.println("The " + i + "th mole's InchiKey is different from the InchiKey generated by our tool");
					System.out.println("Origin title: " + checkInchikey);
					System.out.println("OurGen title: " + inchikey1);
					oldProperties.put("InchiKey",inchikey1);
				}
			}
			if(oldProperties.keySet().contains("SMILES")){
				String checkSmiles = (String) oldProperties.get("SMILES");
				if(!checkSmiles.equals(smile1)){
					System.out.println("The " + i + "th mole's InchiKey is different from the InchiKey generated by our tool");
					System.out.println("Origin title: " + checkSmiles);
					System.out.println("OurGen title: " + smile1);
					oldProperties.put("SMILES",smile1);
				}
			}
						
			CypPre.putAll(oldProperties);
			outMole.addProperties(CypPre);
			//Put it into the resultIAtomContainerSet
			resultMole.addAtomContainer(outMole);
		}
		return resultMole;
	}
	
	/**
	 * Write input molecules with predicted results into the target csv file
	 * @param rawMolecules: The original input compounds
	 * @param outSdfPath: The path to the target csv file
	 * @param predictedResult: predicted results for all compounds
	 * @throws Exception
	 */
	public void outPutCsv(IAtomContainerSet rawMolecules, String outPath, ArrayList<HashMap<String,String>> predictedResult) throws Exception{
		IAtomContainerSet resultMole = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		String csvOutput = outPath;
		FileWriter csvWriter  = new FileWriter(new File(csvOutput));
		//Write the first line: Inchiky, SMILEs, Titile/MoleCounter, PredictedResults
		String titleLine = "Inchiky, SMILES, Title, 1A2, 2B6, 2A6, 2C8, 2C9, 2C19, 2D6, 2|E1, 3A4\n";
		csvWriter.write(titleLine);
		for(int i = 0; i < rawMolecules.getAtomContainerCount(); i++){
			
			IAtomContainer oneMole = rawMolecules.getAtomContainer(i);
			IAtomContainer outMole = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
			outMole = oneMole.clone();
			outMole.setProperties(null);
			//Create proper list
			Map<Object, Object> CypPre = new LinkedHashMap<Object, Object>();
			//Add inchikey
			// Generate factory - throws CDKException if native code does not load
			InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
			// Get InChIGenerator
			InChIGenerator gen1 = factory.getInChIGenerator(oneMole);
			String inchikey1 = gen1.getInchiKey();
					
			//Add SMILEs
			SmilesGenerator sg = new SmilesGenerator().unique();
			String smile1 = sg.create(oneMole);
			
			//Add Title
			String title1orCounter = "";
			String title1 = oneMole.getProperty(CDKConstants.TITLE); 
			if(title1 != null && !title1.isEmpty()){
				title1orCounter =  title1;
			}
			else title1orCounter = "molecule:" + (i+1);
			
			CypPre.put("1A2", predictedResult.get(i).get("1A2"));
			CypPre.put("2B6", predictedResult.get(i).get("2B6"));
			CypPre.put("2A6", predictedResult.get(i).get("2A6"));
			CypPre.put("2C8", predictedResult.get(i).get("2C8"));
			CypPre.put("2C9", predictedResult.get(i).get("2C9"));
			CypPre.put("2C19",predictedResult.get(i).get("2C19"));
			CypPre.put("2D6", predictedResult.get(i).get("2D6"));
			CypPre.put("2E1", predictedResult.get(i).get("2E1"));
			CypPre.put("3A4", predictedResult.get(i).get("3A4"));
			
			if(title1orCounter.contains(",")){
				title1orCounter = "\"" + title1orCounter + "\"";
			}
			String oneSampleLine = inchikey1 + "," + smile1 + "," + title1orCounter + 
					"," + predictedResult.get(i).get("1A2") + "," + predictedResult.get(i).get("2B6") + "," + predictedResult.get(i).get("2A6") + 
					"," + predictedResult.get(i).get("2C8") + "," + predictedResult.get(i).get("2C9") + "," + predictedResult.get(i).get("2C19") +
					"," + predictedResult.get(i).get("2D6") + "," + predictedResult.get(i).get("2E1") + "," + predictedResult.get(i).get("3A4") + "\n"; 
			//Put it into the resultIAtomContainerSet			
			//write the molecule into csv file
			csvWriter.write(oneSampleLine);
			
		}

		csvWriter.close();
	}		
	/**
	 * Reformat the instances of the input molecules so they can be input to the learned models by applying feature selection and normalization
	 * @param testSet: The instances of the user's input molecules
	 * @param attList: The attributes that are used in the learned model
	 * @param meanList: The mean values of each attribute
	 * @param maxList: The max values of each attribute
	 * @param minList: The min values of each attribute
	 * @return
	 * @throws Exception
	 */
	public Instances matchAttributes(Instances testSet, ArrayList<String> attList, ArrayList<String> meanList, ArrayList<String> maxList, ArrayList<String> minList) throws Exception{
		//Create the class atrribute
		List my_nominal_values = new ArrayList(3); 
		my_nominal_values.add("R"); 
		my_nominal_values.add("N"); 
		Attribute lastAttribute  = new Attribute(attList.get(attList.size()-1), my_nominal_values);
		
		ArrayList<Attribute> atts = new ArrayList<Attribute>();
		for(int i = 0; i<attList.size();i++){
			if(i<attList.size()-1){
				Attribute Attribute = new Attribute(attList.get(i));
				atts.add(Attribute);
			}
			else{
				
				atts.add(lastAttribute);
			}
		}
		
		int length = testSet.size();
		int numAttributes = attList.size();
		Instances matched = new Instances("Rel", atts, length);
		Instance sample = new DenseInstance(numAttributes);
		
		for(int j = 0; j< length; j++){
			ArrayList misFeatures = new ArrayList();
			Instance temp = testSet.get(j);
			for(int vidx = 0; vidx < numAttributes; vidx++){
				  Attribute att = atts.get(vidx);
				  if(vidx<numAttributes-1){
					  
					  String stAtt = att.toString();
					  String[] whatever = stAtt.split(" ");
					  
					  if(whatever.length!=3){
						  String combineSpaces = "";
						  for(int kk = 1; kk<whatever.length-1; kk++){
							  if(kk<whatever.length - 2){
								  combineSpaces = combineSpaces + whatever[kk]+ " ";
							  }
							  else combineSpaces = combineSpaces + whatever[kk];
						  }
						  whatever[1] = combineSpaces;
					  }
					 
					  if(whatever[1].contains("'")){
						  whatever[1] = whatever[1].replace("'", "");
					  }
					  Attribute convAtt = testSet.attribute(whatever[1]);
					  double vle = temp.value(convAtt);
					  if(!Double.isNaN(vle)){
						  double min = Double.parseDouble(minList.get(vidx));
						  double max = Double.parseDouble(maxList.get(vidx));
						  double normedVal = (vle-min)/(max-min);
						  sample.setValue(att, normedVal);
					  }
					  else{
						  //If the value of a attribute is not a number, use the mean value of that attribute to replace it.
						  misFeatures.add(whatever[1]);
						  double syn = Double.parseDouble(meanList.get(vidx));
						  sample.setValue(att, syn);

					  }
				  }
				  else{
					  //Set class with random value. Does not matter
					  sample.setValue(lastAttribute, "N");
					  
				  }
			}
			//Output those missing features
			if(misFeatures.size()>0){
				int numMisFeatures = misFeatures.size();
				String misFeatureString = "Molecule" + j + ": " +  numMisFeatures + 
						" features could not be computed from the structure representation, and were thus synthesized: " + misFeatures.get(0);
				for(int i = 1; i < numMisFeatures; i++){
					misFeatureString = misFeatureString + "," + misFeatures.get(i);
				}
				//System.out.println(misFeatureString);
			}
			
			
			matched.add(sample);
		}
		
		return matched;
		
	}
	/**
	 * Initialize the  predictedResults arraylist
	 * @param predictedResult
	 * @param numOfMoles
	 * @return
	 */
	public ArrayList<HashMap<String,String>> initPreResults(ArrayList<HashMap<String,String>> predictedResult, int numOfMoles){
		ArrayList<HashMap<String,String>> afterInit = new ArrayList<HashMap<String,String>>();
		for(int i = 0; i < numOfMoles; i++){
			HashMap<String,String> hm = new HashMap<String,String>();
			hm.put("1A2", "null");
			hm.put("2A6", "null");
			hm.put("2B6", "null");
			hm.put("2C8", "null");
			hm.put("2C9", "null");
			hm.put("2C19", "null");
			hm.put("2D6", "null");
			hm.put("2E1", "null");
			hm.put("3A4", "null");
			afterInit.add(hm);
			
		}
		
		return afterInit;
	}

}

