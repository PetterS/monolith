import unittest

import minimum.core.immutable as immutable


class TestFreeze(unittest.TestCase):
	def test_1d(self):
		assert immutable.freeze([1, 2, 3]) == (1, 2, 3)

	def test_2d(self):
		assert immutable.freeze([1, [2, 3], "4", [[5]]]) == (1, (2, 3), "4",
		                                                     ((5, ), ))

	def test_set(self):
		s = set([1, 2, 3])
		fs = immutable.freeze(set([1, 2, 3]))
		assert s == fs
		assert isinstance(fs, frozenset)

	def test_tuple(self):
		assert immutable.freeze((1, 1)) == (1, 1)
		assert immutable.freeze((1, [5])) == (1, (5, ))

	def test_basic_types(self):
		assert immutable.freeze(1 == 1) == True
		assert immutable.freeze(1.75) == 1.75
		assert immutable.freeze("Petter") == "Petter"

	def test_simple_dict(self):
		d = {1: "a", 2: "b"}
		fd = immutable.freeze(d)
		assert 1 in fd
		assert fd[2] == "b"
		assert len(fd) == 2
		assert immutable.freeze(fd) is fd
		assert hash(fd) == hash(immutable.freeze(d))
		assert fd == fd
		assert fd == immutable.freeze(d)
		with self.assertRaises(TypeError):
			fd.petter = 2
		with self.assertRaises(TypeError):
			fd[3] = "c"

		d2 = {}
		for k, v in fd.items():
			d2[k] = v
		assert immutable.freeze(d2) == fd

		keys = []
		for k in fd:
			keys.append(k)
		assert set(keys) == d.keys()
		assert set(fd.keys()) == set(d.keys())
		assert set(fd.values()) == set(d.values())

	def test_complex_dict(self):
		d = {1: "a", 2: {21: "c", 22: ["d", "e"]}}
		fd = immutable.freeze(d)
		assert type(fd) == type(fd[2])
		assert fd[2][21] == "c"
		assert fd[2][22] == ("d", "e")
		assert len(fd) == 2
		assert immutable.freeze(fd) is fd
		assert hash(fd) == hash(immutable.freeze(d))
		assert fd == fd
		assert fd == immutable.freeze(d)

	def test_dict_str(self):
		assert str(immutable.freeze({})) == "{}"
		assert repr(immutable.freeze({})) == "_FrozenDict({})"
		assert "{0}OK".format(immutable.freeze({})) == "{}OK"

	def test_dict_compare(self):
		d = {"a": 2}
		fd = immutable.freeze(d)
		assert fd == fd
		assert fd == d
		assert d == fd
		assert not fd != fd
		assert not fd != d
		assert not d != fd

	def test_dict_del(self):
		fd = immutable.freeze({"a": 2})
		with self.assertRaises(TypeError):
			del fd._value

	def test_dict_reassign_private(self):
		fd = immutable.freeze({"a": 2})
		with self.assertRaises(TypeError):
			fd._value = {}

	def test_dict_mutate_private(self):
		fd = immutable.freeze({"a": 2})
		with self.assertRaises(TypeError):
			fd._value


if __name__ == "__main__":
	unittest.main()
