# Getting Started

This guide will help you quickly install `slow-rotations`. A molecular dynamics simulation analysis package to flag slow torsional motions that may indicate sampling issues. 

## Installation

1. Clone the [github repository](https://github.com/MobleyLab/slow-rotations/tree/main#)
	```
	git clone git@github.com:MobleyLab/slow-rotations.git
	```

2. Navigate to the top level of the slow-rotations github repository

2. Install the environment needed to run slow-rotations
	```
	mamba env install -f environment.yaml
	```

3. Install slow-rotations from the top level of the slow-rotations github repository
	```
	pip install -e .
	```

3. You should be able to run `slow-rotations` scripts now!


## Run an example for small molecules
Example scripts and outputs can be found in the `slow-rotations/examples` path in the github directory

1. Download the example trajectories 

2. Change the file paths to trajectories in `slow-rotations/examples/example_smallmolecule.py` to match the path where you have stored the example trajectories

3. Run the small molecule example
	```
	python examples/example_smallmolecule.py
	```

4. Look for the output torsion images and `.json` file in `slow-rotations/examples`

	![Example](./images/8_6_3_0.png)

5. Analyze the each torsion for potential sampling issues
	```
	python examples/example_analysis_smallmolecule.py
	```

6. Check `examples/mol_flaggedtorsions.csv` to see which torsions are flagged.
	* In the example below all 3 repeats are flagged for low transitions and repeat 1 is flagged for a missing state (see step 4)

	| residue | torsion       | torsion_sys               | repeat | low transitions | missing states |
	|---------|---------------|---------------------------|--------|----------------|----------------|
	| LIG1    | [8, 6, 3, 0]  | [2554, 2552, 2549, 2546] | 0      | True           |                |
	| LIG1    | [8, 6, 3, 0]  | [2554, 2552, 2549, 2546] | 1      | True           | 1              |
	| LIG1    | [8, 6, 3, 0]  | [2554, 2552, 2549, 2546] | 2      | True           |                |

## Run an example for protein sidechains
Example scripts and outputs can be found in the `slow-rotations/examples` path in the github directory

1. Download the example trajectories 

2. Change the file paths to trajectories in `slow-rotations/examples/example_sidechain.py` to match the path where you have stored the example trajectories

3. Run the small molecule example
	```
	python examples/example_sidechain.py
	```

4. Look for the output torsion images and `.json` file in `slow-rotations/examples`

	![Example](./images/2291_2292_2295_2296.png)

5. Analyze the each torsion for potential sampling issues
	```
	python examples/example_analysis_sidechain.py
	```

6. Check `examples/sidechain_flaggedtorsions.csv` to see which torsions are flagged.
	* In the example below all 3 repeats are flagged for low transitions and repeat 1 is flagged for a missing state (see step 4)

	| residue | chi | torsion               | repeat | low transitions | missing states |
	|---------|-----|----------------------|--------|----------------|----------------|
	| VAL154  | 1   | [2291, 2292, 2295, 2296] | 0      | True           |                |
	| VAL154  | 1   | [2291, 2292, 2295, 2296] | 1      | True           | 0              |







