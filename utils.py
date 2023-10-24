import warnings

def warn(msg: str):
	warnings.warn(msg)


class NoMoleculeError(Exception):
	pass

class NotImplementedError(Exception):
	pass

