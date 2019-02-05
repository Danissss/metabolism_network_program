package xuan.drug_porter.descriptorUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Point3d;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.qsar.DescriptorEngine;
import org.openscience.cdk.qsar.IAtomicDescriptor;
import org.openscience.cdk.qsar.IDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.IPAtomicLearningDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor;
import org.openscience.cdk.qsar.result.BooleanResult;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.qsar.result.IntegerArrayResult;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.openscience.cdk.tools.HOSECodeGenerator;

public class GetAtomicDescriptors {
	

	static List<String> classNames = DescriptorEngine.getDescriptorClassNameByPackage("org.openscience.cdk.qsar.descriptors.atomic",
            null);
	private static  DescriptorEngine ENGINE = new DescriptorEngine(classNames, null);

	//	private static DescriptorEngine ENGINE = new DescriptorEngine(DescriptorEngine.ATOMIC);
	// Find nearest atom to all atoms in a molecule
	// 
	public static ArrayList<ArrayList<String>> getNearestAtoms(IAtomContainer mole) {
		
		ArrayList<ArrayList<String>> atom_distances = new ArrayList<ArrayList<String>>();

		int atomCount = mole.getAtomCount();

		List<IAtom> atoms = new ArrayList<IAtom>();
		for (int i = 0; i < atomCount; i++) {
			atoms.add(mole.getAtom(i));
		}

		for (int i = 0; i < atoms.size(); i++) {
			Double[] distances = new Double[atoms.size()];
				
			for (int j = 0; j < atoms.size(); j++) {
				if (j == i) {
						// Large number so that sorting puts it last
					distances[j] = 99999.12;
					continue;
				}
				
				Point3d firstPoint = atoms.get(i).getPoint3d();
				Point3d secondPoint = atoms.get(j).getPoint3d();
				Double distance = firstPoint.distance(secondPoint);
				distances[j] = distance;

			}

			// put the nearest atom at front
			ArrayList<String> indices = new ArrayList<String>();
			Double[] d = distances.clone(); // clone the original list (unsorted)
			Arrays.sort(d); // sort;
			List<Double> d_list = Arrays.asList(distances); // put into

			for (int j = 0; j < distances.length; j++) {
				String index = String.valueOf(d_list.indexOf(d[j])); // get index of that changed atoms
				indices.add(index);
			}

//			it return the nearest atom, not the actual distance
			atom_distances.add(indices);
				
			
			}
		return atom_distances;
	}
		


	/**
	 * Calculate descriptors. Omits IPMolecularLearningDescriptor
	 *
	 * @param string
	 *            path to SDF input file
	 * @param string
	 *            path to CSV output file
	 * @param string
	 *            comma-seperated list of descriptor names (if empty, all
	 *            descriptors will be calculated)
	 */
	/*
	 *  this is used for the NmrExperiment
	 *  get all atomic descriptor for NmrExperiment
	 */
	public static ArrayList<Double[]> getAtomicDescriptor(IAtomContainer mole, List<IAtom> atoms,String descNamesStr) throws java.io.IOException {
		
	
		List<IDescriptor> descriptors = ENGINE.getDescriptorInstances();
		List<String> descNames = Arrays.asList(descNamesStr.split(","));
		ArrayList<String> colNames = new ArrayList<String>();
		ArrayList<Double[]> values = new ArrayList<Double[]>();
		
		
		// get each descriptors
//		int num_descriptor = 0;
		for (IDescriptor desc : descriptors) {
			if (desc instanceof IPAtomicLearningDescriptor)
				continue;
			String tname = desc.getClass().getName();
			String[] tnamebits = tname.split("\\.");
			tname = tnamebits[tnamebits.length - 1];
			if ((descNamesStr.length() > 0) && (!descNames.contains(tname)))
				continue;
			String[] colNamesArr = desc.getDescriptorNames();
			for (int idx = 0; idx < colNamesArr.length; idx++) {
				colNamesArr[idx] = tname + "-" + colNamesArr[idx];
			}

			colNames.addAll(Arrays.asList(colNamesArr));
			
			
				
				
			
			//doesn't need to calculate all the descriptor:
			
			
			// call the computeListsAtomic to get the desired atomic value 
//			int atomCount = mole.getAtomCount();
//			List<IAtom> atoms = new ArrayList<IAtom>();
//			for (int i = 0; i < atomCount; i++) {
//				atoms.add(mole.getAtom(i));
//			}
			
			// each row is the descriptor value of desc for each atoms; i.e. if there are 4 atom, then, each row has four value;
			// totally number of descriptor rows.
			values.addAll(computeDescriptorsAtomic(mole, atoms, (IAtomicDescriptor) desc));
			try {
				getHoseCodesForMolecule(mole); // return: ArrayList<String>
			}
			catch (Exception e){
				System.out.println(mole.getTitle());
			}
//			num_descriptor++;
				
		}
//		System.out.println("Number of Descriptors: " + num_descriptor); // 29
		return values;
	}

	/*
	 * for computeListAtomic 
	 * Most important function (thank god, finally)
	 * the value for descriptor comes from this function call
	 * input: atomContainer mol; list of atoms; descriptors
	 * calculate each descriptor for each atoms (23 atoms)
	 */
	public static List<Double[]> computeDescriptorsAtomic(IAtomContainer mol, List<IAtom> atoms,
			IAtomicDescriptor descriptor) {
		List<Double[]> vv = new ArrayList<Double[]>();

//		// total 23 atoms (same as sdf file shows)

		vv.add(new Double[atoms.size()]);
		
		
		// iterate each atom
		for (int i = 0; i < atoms.size(); i++) {
			if (atoms.get(i) == null) {
				vv.get(0)[i] = null;
			} else {
				try {
					IDescriptorResult res = descriptor.calculate(atoms.get(i), mol).getValue();
					//System.out.println(res.toString());// res contain all the value for each atom
					if (res instanceof IntegerResult) {
						vv.get(0)[i] = (double) ((IntegerResult) res).intValue();
//						 System.out.println("IntegerResult"+vv.get(0)[i]);
					} else if (res instanceof DoubleResult) {
						vv.get(0)[i] = ((DoubleResult) res).doubleValue();
//						 System.out.println("DoubleResult"+vv.get(0)[i]);
					} else if (res instanceof DoubleArrayResult) {
						vv.get(0)[i] = ((DoubleArrayResult) res).get(0);
//						 System.out.println("DoubleArrayResult"+vv.get(0)[i]);
					} else if (res instanceof IntegerArrayResult) {
						vv.get(0)[i] = (double) ((IntegerArrayResult) res).get(0);
//						 System.out.println("IntegerArrayResult"+vv.get(0)[i]);
					} else if (res instanceof BooleanResult) {
						String result_bool = res.toString();
						if (result_bool == "false") {
							vv.get(0)[i] = 0.0;
						} else {
							vv.get(0)[i] = 1.0;
						}
					}else{
						throw new IllegalStateException(
								"Unknown idescriptor result value for '" + descriptor + "' : " + res.getClass());
					}
				}catch (Throwable e) {
					System.err.println("Could not compute cdk feature " + descriptor);
					e.printStackTrace();
					vv.get(0)[i] = 0.0;
				}
			}
//			System.out.println(atoms.get(i).getSymbol());

//			System.out.println(Arrays.toString(vv.get(0)));
			
			if (vv.get(0)[i] != null && (vv.get(0)[i].isNaN() || vv.get(0)[i].isInfinite()))
				vv.get(0)[i] = 0.0;
		}
		
		
		return vv;
	}
	
	
	public static ArrayList<String> getHoseCodesForMolecule(IAtomContainer mol) {
		HOSECodeGenerator hoseG = new HOSECodeGenerator();
		ArrayList<String> hoseCodes = new ArrayList<String>();

		int atomCount = mol.getAtomCount();
		List<IAtom> atoms = new ArrayList<IAtom>();
		for (int i = 0; i < atomCount; i++) {
			try {
				String hose = hoseG.getHOSECode(mol, mol.getAtom(i), 0);
				hoseCodes.add(hose);
//				System.out.println("HOSE = " + hose + "\n");
			} catch (CDKException e) {
				e.printStackTrace();
			}
		}
		return hoseCodes;
	}
}
