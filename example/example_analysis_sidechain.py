from  slow_rotations.torsiondata import *
from  slow_rotations.analysis import *
import csv


torsion_data_ifile = "sidechain_torsiondata.json"
flagged_torsions_ofile = "sidechain_flagged_torsions.csv"


with open(torsion_data_ifile, "r") as f:
    json_str = f.read()

data = TorsionData.from_json(json_str)

problem_torsions = []

for tname in data.list_torsions():
	torsion = data.get_torsion(tname)

	for rnum,rpt in torsion.repeats.items():
		repeat_info = {
			'residue': rpt.residue,
			'chi': rpt.chi,
			'torsion':  rpt.torsion_idx,
			'repeat': rnum, 
			'low transitions': False,
			'missing states': '',
		}
		found_problem = False
		if not rpt.symmetry:
			low_transitions = check_transitions(rpt)
			if check_transitions(rpt):
				repeat_info['low transitions'] = True
				found_problem = True

			missing = check_states(rpt)
			if len(missing) > 0:
				repeat_info['missing states'] = " ".join(map(str, missing))
				found_problem = True

		if found_problem:
			problem_torsions.append(repeat_info)


with open(flagged_torsions_ofile, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=problem_torsions[0].keys())
    writer.writeheader()
    writer.writerows(problem_torsions)
