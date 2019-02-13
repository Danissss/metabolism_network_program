package xuan.drug_porter.reactant;

import java.util.HashMap;
import java.util.Map;

import weka.classifiers.Classifier;
import weka.core.Instances;

public class phase2react {

	
	
	protected static String current_dir = System.getProperty("user.dir");
	private static HashMap<String, String> substrate_model_path = create_substrate_model_path("substrate");
	
	
	
	
	/**
	 * return all transporter model result
	 * @param instance
	 * @return
	 * @throws Exception 
	 */
	
	public static HashMap<String,String> RunClassifierForAll(Instances instance) throws Exception {
		
		HashMap<String,String> result =  new HashMap<String,String>();
		String[] model_name = new String[] {"UGT","SULT","NAT","GST","COMT"};
		
		for(int i=0 ; i<model_name.length; i++) {
			String s_model_name = model_name[i];
			HashMap<String,String> getback_result = RunClassifier(instance,s_model_name);
			for (Map.Entry<String,String> entry : getback_result.entrySet()){
				result.put(entry.getKey(), entry.getValue());
			}
		}
		
		
		return result;
		
		
	}
	
	/**
	 * Get single transporter model 
	 * @param instance
	 * @param enzyme_name
	 * @return
	 * @throws Exception
	 */
	public static HashMap<String,String> RunClassifier(Instances instance, String enzyme_name) throws Exception{
		
		HashMap<String,String> classified_result = new HashMap<String,String>();
		Classifier model = (Classifier) weka.core.SerializationHelper.read(substrate_model_path.get(enzyme_name));
		double result = model.classifyInstance(instance.get(0));		
		String key_string = String.format("%s_substrate", enzyme_name);
		
		if(result == 0.0) {
			classified_result.put(key_string, "non-substrate");
		}else {
			classified_result.put(key_string, "substrate");
		}
		
		
		
		return classified_result;
	}
	
	
	/**
	 * 
	 * @param model_type
	 * @return
	 */
	public static HashMap<String, String> create_substrate_model_path(String model_type){
		
		HashMap<String,String> model_map = new HashMap<String,String>();
		model_map.put("UGT", String.format("%s/Model/PhaseIIReactModel/%s_RANDOMFOREST_SUBSTRATE.model", current_dir,"UGT"));
		model_map.put("SULT", String.format("%s/Model/PhaseIIReactModel/%s_RANDOMFOREST_SUBSTRATE.model", current_dir,"SULT"));
		model_map.put("NAT", String.format("%s/Model/PhaseIIReactModel/%s_RANDOMFOREST_SUBSTRATE.model", current_dir,"NAT"));
		model_map.put("GST", String.format("%s/Model/PhaseIIReactModel/%s_RANDOMFOREST_SUBSTRATE.model", current_dir,"GST"));
		model_map.put("COMT", String.format("%s/Model/PhaseIIReactModel/%s_RANDOMFOREST_SUBSTRATE.model", current_dir,"COMT"));
		return model_map;
	}
	
}
