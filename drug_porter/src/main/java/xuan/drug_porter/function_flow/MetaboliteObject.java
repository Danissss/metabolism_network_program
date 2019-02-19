package xuan.drug_porter.function_flow;

import java.util.ArrayList;

//import org.openscience.cdk.interfaces.IAtomContainer;

public class MetaboliteObject {
	
	
	public String SMILES = new String();
	public String ParentSMILES = new String();
	public ArrayList<String> ChildSMILES = new ArrayList<String>();
	public int NumberOfChild = new Integer(0);
	public ArrayList<String> TargetSubstrate = new ArrayList<String>();
	public ArrayList<String> TargetInhibitor = new ArrayList<String>();
	public boolean bioactive = false;
 	
	
	
	
	/**
	 * constructor 
	 * @param smiles
	 */
	public MetaboliteObject(String smiles) {
		this.SMILES = smiles;
		
	}
	
	
	/**
	 * constructor 
	 * @param smiles
	 * @param parentsSmiles
	 */
	public MetaboliteObject(String smiles, String parentsSmiles) {
		this.SMILES = smiles;
		this.ParentSMILES = parentsSmiles;
		
		
	
	}
	
	/**
	 * constructor
	 * @param smiles
	 * @param parentsSmiles
	 * @param child
	 */
	public MetaboliteObject(String smiles, String parentsSmiles, ArrayList<String> child) {
		this.SMILES = smiles;
		this.ParentSMILES = parentsSmiles;
		this.ChildSMILES = child;
		if(this.ChildSMILES.size()!=0) {
			UpdateNumberOfChild(this.ChildSMILES.size());
		}
		
	}
	
	
	/**
	 * Set BioActive as true
	 */
	public void BioActive(boolean what) {
		this.bioactive = what;
	}
	
	/**
	 * each molecule can only have one parents
	 * @param container
	 * @param smiles
	 */
	public void AddParent(String smiles) {
		this.ParentSMILES = smiles;
		
	}
	
	/**
	 * add child smiles 
	 * @param smiles
	 */
	public void AddChild(String smiles) {
		this.ChildSMILES.add(smiles);
		UpdateNumberOfChild(this.ChildSMILES.size());
		
	}
	
	/**
	 * get number of child 
	 * @return
	 */
	public Integer GetNumOfChild() {
		return this.ChildSMILES.size();
	}
	
	/**
	 * return child in arraylist
	 * @return
	 */
	public ArrayList<String> GetChild(){
		return this.ChildSMILES;
	}
	
	/**
	 * return parents smiles 
	 * @return
	 */
	public String GetParent() {
		return this.ParentSMILES;
	}
	
	/**
	 * 
	 * @return
	 */
	public ArrayList<String> GetAllTargetSubstrate(){
		return TargetSubstrate;
	}
	
	/**
	 * 
	 * @return
	 */
	public ArrayList<String> GetAllTargetInhibitor(){
		return TargetInhibitor;
	}
	
	
	/**
	 * 
	 * @return
	 */
	public void AddTargetAsSubstrate(String TargetName){
		
		if(!this.TargetSubstrate.contains(TargetName)) {
			this.TargetSubstrate.add(TargetName);
		}
		
	}
	
	
	/**
	 * 
	 * @return
	 */
	public void AddTargetAsInhibitor(String TargetName){
		
		if(!this.TargetInhibitor.contains(TargetName)) {
			this.TargetInhibitor.add(TargetName);
		}
	}
	
	/**
	 * 
	 * @param TargetName
	 */
	public void RemoveTargetFromSubstrate(String TargetName) {
		if (this.TargetSubstrate.contains(TargetName)){
			this.TargetSubstrate.remove(TargetName);
		}
	}
	
	
	/**
	 * 
	 * @param TargetName
	 */
	public void RemoveTargetFromInhibitor(String TargetName) {
		if (this.TargetInhibitor.contains(TargetName)){
			this.TargetInhibitor.remove(TargetName);
		}
	}
	
	
	
	/**
	 * 
	 * @param number
	 */
	private void UpdateNumberOfChild(int number) {
		this.NumberOfChild = number;
		
	}
	
	/**
	 * 
	 * @return
	 */
	boolean IsRoot() {
		if (this.SMILES == this.ParentSMILES) {
			return true;
		}else {
			return false;
		}
	}
	
	/**
	 * 
	 * @return
	 */
	String GetSMILES() {
		return this.SMILES;
	}
	

}
