package xuan.drug_porter.descriptorUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.LinkedHashMap;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.qsar.DescriptorEngine;
import org.openscience.cdk.qsar.IAtomicDescriptor;
import org.openscience.cdk.qsar.IDescriptor;
import org.openscience.cdk.qsar.IMolecularDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.IPAtomicLearningDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.qsar.result.IntegerArrayResult;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk.fingerprint.KlekotaRothFingerprinter;

/**
 * Hello world!
 *
 */
public class GenerateFeatureSingle 
{
	
	
	List<String> classNames = DescriptorEngine.getDescriptorClassNameByPackage("org.openscience.cdk.qsar.descriptors.molecular",
            null);
	DescriptorEngine descriptoEngine = new DescriptorEngine(classNames, null);
	
	/**
	 * Given an IAtomContainer of a molecule, generate a string that contains all raw feature values for that molecule
	 * @param IAtomContainer molecule
	 * @return IAtomContainerSet that contains all molecules in the sdf file        
	 * @throws Exception
	 * @author Siyan Tian, Xuan
	 */
	
	
	public String generateOneinstance(IAtomContainer mole,String featureType ) throws Exception {
		StringBuffer sb = new StringBuffer();
		ChemSearcher cs = new ChemSearcher();
		PubchemFingerprinter pbf 	= new PubchemFingerprinter(SilentChemObjectBuilder.getInstance());
		MACCSFingerprinter maccs 	=  new MACCSFingerprinter(SilentChemObjectBuilder.getInstance());
		
		LinkedHashMap<String, String> fpatterns = cs.getRINFingerprintPatterns();
		FeatureGeneration fgen = new FeatureGeneration();
		
		IAtomContainer container = mole;
		
		IAtomContainer prepContainer = MoleculeExplorer.preprocessContainer(container);
		String[] gg = fgen.generateExtendedMolecularFeatures(prepContainer).split(",");
		
		//this is molecular featuresm
		String extendedFeatures = StringUtils.join(fgen.generateExtendedMolecularFeatures(prepContainer).split(","), "\t");
		String molecularFeatures = StringUtils.join(fgen.generateExtendedMolecularFeatures(prepContainer).split(","), ",");
		
		// nonBitFeature contains the feature that don't have bit feature
		String nonBitFeature = extendedFeatures;
		String[] nonBitFeatures = nonBitFeature.split("\t");
		
		ArrayList<Double> bioTransformerFingerprint_bits = cs.generateClassyfireFingerprintAsDouble(prepContainer, fpatterns).getBitValues();
		
		//print bioTransformerFingerprint_bits separated by comma
		String bioTFinger_bits = "";
		for(int x = 0; x < bioTransformerFingerprint_bits.size(); x++){
			bioTFinger_bits =  bioTFinger_bits + String.valueOf(bioTransformerFingerprint_bits.get(x)) + ",";
		}
		
		//bioTFinger_bits
		bioTFinger_bits = bioTFinger_bits.substring(0, bioTFinger_bits.length()-1);
		
		
		//extendedFeatures = molecular Features + bioTFinger_bits
		for(int x = 0; x < bioTransformerFingerprint_bits.size(); x++){
			extendedFeatures =  extendedFeatures + "\t" + String.valueOf(bioTransformerFingerprint_bits.get(x));

		}
		
		
		ArrayList<Double> fingerprint_bits = new ArrayList<Double>();
		IBitFingerprint fingerp	= pbf.getBitFingerprint(prepContainer);

		int[] onbits = fingerp.getSetbits();

		for(int kp = 0; kp < 881; kp++){
			fingerprint_bits.add(0.0);
		}
		for(int o = 0; o < onbits.length; o++){
			fingerprint_bits.set(onbits[o], 1.0);
		}
		
		String pubchemFingerPrint = "";
		pubchemFingerPrint = pubchemFingerPrint + "," + StringUtils.join(fingerprint_bits,",");
		
		extendedFeatures =  extendedFeatures + "\t" + StringUtils.join(fingerprint_bits,"\t");
		
			
		ArrayList<Double> maccs_fingerprint_bits = new ArrayList<Double>();
		IBitFingerprint maccs_fingerp		= maccs.getBitFingerprint(prepContainer);
			
		int[] maccs_onbits = maccs_fingerp.getSetbits();
			
		for(int kp = 0; kp < 166; kp++){
			maccs_fingerprint_bits.add(0.0);
		}
		for(int o = 0; o < maccs_onbits.length; o++){
			maccs_fingerprint_bits.set(maccs_onbits[o], 1.0);
		}
		
		//System.out.println("4::"+extendedFeatures);
		String maccFingerprint = "";
		maccFingerprint = maccFingerprint + "," + StringUtils.join(maccs_fingerprint_bits,",");
		extendedFeatures =  extendedFeatures + "\t" + StringUtils.join(maccs_fingerprint_bits,"\t");
		
		String finalFeatureValues = extendedFeatures;
		String[] temp = extendedFeatures.split("\t");
		
		
		//select which to return:
		if(featureType == "fingerprint") {
			//1197
			return bioTFinger_bits;
		}
		else if (featureType == "pubchem"){
			return pubchemFingerPrint;
		}
		else if (featureType == "macc") {
			return maccFingerprint;
		}
		else if (featureType == "ALL") {
			
			String AllFeature = molecularFeatures +","+ bioTFinger_bits+pubchemFingerPrint+ maccFingerprint;
			return AllFeature;
		}
		else {
			return molecularFeatures;
		}
		
	
	}
	/**
	 * 
	 * @param mole
	 * @param featureType useless parameters
	 * @param returnType: value (descriptor value); name (descriptor name);
	 * @return depends on returnType 
	 * @throws Exception
	 * @author xuan
	 */
	public String generate_all_cdk_molecular_descriptor(IAtomContainer mole,String featureType, String returnType) throws Exception {
		
		//List<String> descriptor_name = descriptoEngine.getDescriptorClassNames();
		List<IDescriptor> descriptors = 	descriptoEngine.getDescriptorInstances();

		
		if (returnType == "value") {
			String temp = null;
			
			for (IDescriptor desc : descriptors) {
				try {
					IDescriptorResult res = ((IMolecularDescriptor) desc).calculate(mole).getValue();
//					int temp_int  = ((IntegerResult) res).intValue();
//					System.out.println(temp_int);
//					System.out.println(res.toString());// res contain all the value for each atom
					
					if (res instanceof IntegerResult) {
						temp = temp+","+res.toString();
						//System.out.println("Result:"+temp);
						//System.out.println("IntegerResult"+String.valueOf(temp));
//						 System.out.println("IntegerResult"+vv.get(0)[i]);
					} else if (res instanceof DoubleResult) {
						temp = temp+","+res.toString();
						//System.out.println("Result:"+temp);
						//System.out.println("IntegerResult"+String.valueOf(temp));
//						 System.out.println("DoubleResult"+vv.get(0)[i]);
					} else if (res instanceof DoubleArrayResult) {
						temp = temp+","+res.toString();
						//System.out.println("Result:"+temp);
						//System.out.println("IntegerResult"+String.valueOf(temp));
//						 System.out.println("DoubleArrayResult"+vv.get(0)[i]);
					} else if (res instanceof IntegerArrayResult) {
						temp = temp+","+res.toString();
						//System.out.println("Result:"+temp);
						//System.out.println("IntegerResult"+String.valueOf(temp));
//						 System.out.println("IntegerArrayResult"+vv.get(0)[i]);
					} else
						throw new IllegalStateException(
								"Unknown idescriptor result value for '" + desc + "' : " + res.getClass());
				} catch (Throwable e) {
					System.err.println("Could not compute cdk feature " + desc);

				}
				
			}
//			String[] attribute_name_list = temp.split(",");
//			System.out.println(attribute_name_list.length);
			return temp;
		}
		else if (returnType == "name") {
			// only print out the feature name
			String attribute_name = "";
			for (IDescriptor desc : descriptors) {
				String[] desr_name = desc.getDescriptorNames();
				String str = String.join(",", desr_name);
				attribute_name = attribute_name + str + ",";
			}

			return attribute_name;
		}
		else {
			System.out.println("Unknown return type");
			System.exit(1);
		}
		
		return null;
		
	}
	
	
}
