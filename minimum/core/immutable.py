from typing import Any


class _FrozenDict:
	__slots__ = "_value", "_hash"

	def __init__(self, value):
		d = {}
		object.__setattr__(self, "_value", d)
		for key, value in value.items():
			d[key] = freeze(value)
		object.__setattr__(self, "_hash", None)

	def items(self):
		return object.__getattribute__(self, "_value").items()

	def keys(self):
		return object.__getattribute__(self, "_value").keys()

	def values(self):
		return object.__getattribute__(self, "_value").values()

	def __iter__(self):
		return iter(object.__getattribute__(self, "_value"))

	def __len__(self):
		return len(object.__getattribute__(self, "_value"))

	def __bool__(self):
		return bool(object.__getattribute__(self, "_value"))

	def __getitem__(self, key):
		return object.__getattribute__(self, "_value")[key]

	def __contains__(self, key):
		return key in object.__getattribute__(self, "_value")

	def __eq__(self, rhs):
		if type(rhs) == dict:
			return object.__getattribute__(self, "_value") == rhs
		elif type(rhs) == _FrozenDict:
			return (object.__getattribute__(
			    self, "_value") == object.__getattribute__(rhs, "_value"))
		else:
			return False

	def __hash__(self):
		existinghash = object.__getattribute__(self, "_hash")
		if not existinghash:
			value = object.__getattribute__(self, "_value")
			newhash = hash(frozenset(value.items()))
			object.__setattr__(self, "_hash", newhash)
			return newhash
		else:
			return existinghash

	def __str__(self):
		return str(object.__getattribute__(self, "_value"))

	def __repr__(self):
		return "_FrozenDict({})".format(
		    repr(object.__getattribute__(self, "_value")))

	def __format__(self, s):
		return object.__getattribute__(self, "_value").__format__(s)

	def __delattr__(self, _):
		raise TypeError("Attribute is immutable")

	def __setattr__(self, name, value):
		raise TypeError("Attribute is inaccessible")

	def __getattribute__(self, name):
		if name != "__slots__" and name in self.__slots__:
			raise TypeError("Attribute is inaccessible")
		return object.__getattribute__(self, name)


def freeze(obj: Any) -> Any:
	"""Recursively makes a data structure immutable.

	Lists become tuples and sets become frozensets. Basic immutable data types
	are kept as-is. Dicts are replaces with a custom frozen class.
	"""
	if type(obj) in (bool, bytes, complex, float, int, range, str,
	                 _FrozenDict):
		return obj
	elif isinstance(obj, (list, tuple)):
		return tuple(freeze(o) for o in obj)
	elif isinstance(obj, set):
		return frozenset(freeze(o) for o in obj)
	elif isinstance(obj, dict):
		return _FrozenDict(obj)
	else:
		# There are certain ways of modifying objects to
		# make them more or less immutable. But there are
		# no guarantees that the object will continue
		# working afterwards.
		raise ValueError("Don't know how to freeze this object.")
