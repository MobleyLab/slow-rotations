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


## Run an example
Example scripts and outputs can be found in the `slow-rotations/examples` path in the github directory

1. Download the example trajectories 

2. Change the file paths to trajectories in `slow-rotations/examples/example_smallmolecule.py` to match the path where you have stored the example trajectories

3. Run the small molecule example
	```
	python examples/example_smallmolecule.py
	```

4. Look for the output torsion images and `.json` file in `slow-rotations/examples`

	![Example](./images/14_7_4_1.png)

5. Analyze the each torsion for potential sampling issues
	```
	python examples/example_analysis.py
	```







