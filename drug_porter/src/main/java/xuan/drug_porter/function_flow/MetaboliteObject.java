package xuan.drug_porter.function_flow;

import java.util.ArrayList;

//import org.openscience.cdk.interfaces.IAtomContainer;

public class MetaboliteObject {
	
	
	public String SMILES = new String();
	public String ParentSMILES = new String();
	public ArrayList<String> ChildSMILES = new ArrayList<String>();
	public int NumberOfChild = new Integer(0);
	
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
	 * @param number
	 */
	private void UpdateNumberOfChild(int number) {
		this.NumberOfChild = number;
		
	}
	

}
