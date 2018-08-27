import io
import struct


def write(f: io.BufferedIOBase, data: bytes) -> None:
	if len(data) < 127:
		f.write(struct.pack('=b', len(data)))
	else:
		f.write(struct.pack('=bq', -1, len(data)))
	f.write(data)


def read(f: io.BufferedIOBase) -> bytes:
	size, = struct.unpack('=b', f.read(1))
	if size < 0:
		size, = struct.unpack('=q', f.read(8))
	return f.read(size)
