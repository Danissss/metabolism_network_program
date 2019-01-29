package xuan.drug_porter.descriptorUtils;



import weka.core.Instances;
import weka.core.converters.ArffSaver;
import weka.core.converters.CSVLoader;
 
import java.io.File;
 
public class ConvertTOArff {
  /**
   * takes 2 arguments:
   * - CSV input file
   * - ARFF output file
   */
	public static void CSVToArff(String filePath,String arff_path) throws Exception {
		
		// load CSV
		CSVLoader loader = new CSVLoader();
		loader.setSource(new File(filePath));
		Instances data = loader.getDataSet();
 
		// save ARFF
		ArffSaver saver = new ArffSaver();
		saver.setInstances(data);
		saver.setFile(new File(arff_path));
		//saver.setDestination(new File(arffFilePath));
		saver.writeBatch();
	}
	
}