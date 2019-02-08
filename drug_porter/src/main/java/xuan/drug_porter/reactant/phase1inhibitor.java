package xuan.drug_porter.reactant;

import java.util.HashMap;
import java.util.Map;

import weka.classifiers.Classifier;
import weka.core.Instances;


public class phase1inhibitor {

	
	
	protected static String current_dir = System.getProperty("user.dir");
	private static HashMap<String, String> inhibitor_model_path = create_inhibitor_model_path("substrate");
	
	
	
	
	/**
	 * return all transporter model result
	 * @param instance
	 * @return
	 * @throws Exception 
	 */
	
	public static HashMap<String,String> run_classifier_for_all_transporter(Instances instance) throws Exception {
		
		HashMap<String,String> result =  new HashMap<String,String>();
		String[] model_name = new String[] {"CYP1A2","CYP2A6","CYP2B6","CYP2C8","CYP2C9","CYP2C19","CYP2C19","CYP2D6","CYP2E1","CYP3A4"};
		
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
	public static HashMap<String,String> run_classifier(Instances instance, String enzyme_name) throws Exception{
		
		HashMap<String,String> classified_result = new HashMap<String,String>();
		Classifier model = (Classifier) weka.core.SerializationHelper.read(inhibitor_model_path.get(enzyme_name));
		double result = model.classifyInstance(instance.get(0));		
		String key_string = String.format("%s_inhibitor", enzyme_name);
		
		if(result == 0.0) {
			classified_result.put(key_string, "non-inhibitor");
		}else {
			classified_result.put(key_string, "inhibitor");
		}
		
		
		
		return classified_result;
	}
	
	
	/**
	 * 
	 * @param model_type
	 * @return
	 */
	public static HashMap<String, String> create_inhibitor_model_path(String model_type){
		
		HashMap<String,String> model_map = new HashMap<String,String>();
		model_map.put("CYP1A2", String.format("%s/Model/PhaseIinhibitorModel/%s_RANDOMFOREST_INHIBITOR.model", current_dir,"CYP1A2"));
		model_map.put("CYP2A6", String.format("%s/Model/PhaseIinhibitorModel/%s_RANDOMFOREST_INHIBITOR.model", current_dir,"CYP2A6"));
		model_map.put("CYP2B6", String.format("%s/Model/PhaseIinhibitorModel/%s_RANDOMFOREST_INHIBITOR.model", current_dir,"CYP2B6"));
		model_map.put("CYP2C8", String.format("%s/Model/PhaseIinhibitorModel/%s_RANDOMFOREST_INHIBITOR.model", current_dir,"CYP2C8"));
		model_map.put("CYP2C9", String.format("%s/Model/PhaseIinhibitorModel/%s_RANDOMFOREST_INHIBITOR.model", current_dir,"CYP2C9"));
		model_map.put("CYP2C19", String.format("%s/Model/PhaseIinhibitorModel/%s_RANDOMFOREST_INHIBITOR.model", current_dir,"CYP2C19"));
		model_map.put("CYP2D6", String.format("%s/Model/PhaseIinhibitorModel/%s_RANDOMFOREST_INHIBITOR.model", current_dir,"CYP2D6"));
		model_map.put("CYP2E1", String.format("%s/Model/PhaseIinhibitorModel/%s_RANDOMFOREST_INHIBITOR.model", current_dir,"CYP2E1"));
		model_map.put("CYP3A4", String.format("%s/Model/PhaseIinhibitorModel/%s_RANDOMFOREST_INHIBITOR.model", current_dir,"CYP3A4"));
		return model_map;
	}
	
}
