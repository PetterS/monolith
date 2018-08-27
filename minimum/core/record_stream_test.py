import io
import struct
import unittest

from minimum.core import record_stream


class RecordStreamTest(unittest.TestCase):
	def test_basic(self):
		empty = b""
		small = b"Petter"
		big = 500 * b"P"

		stream = io.BytesIO(b'')
		record_stream.write(stream, small)
		record_stream.write(stream, small)
		record_stream.write(stream, big)
		record_stream.write(stream, empty)
		record_stream.write(stream, small)
		record_stream.write(stream, big)

		stream.seek(0)

		self.assertEqual(record_stream.read(stream), small)
		self.assertEqual(record_stream.read(stream), small)
		self.assertEqual(record_stream.read(stream), big)
		self.assertEqual(record_stream.read(stream), empty)
		self.assertEqual(record_stream.read(stream), small)
		self.assertEqual(record_stream.read(stream), big)
		self.assertFalse(stream.read(1))

	def test_raise_on_read_failure(self):
		stream = io.BytesIO(b'')
		with self.assertRaises(struct.error):
			record_stream.read(stream)
