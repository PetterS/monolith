import minimum.core

ip_pb2 = minimum.core.load_module("minimum.linear.ip_pb2",
                                  "minimum/linear/ip_pb2.py")
Bound = ip_pb2.Bound
Constraint = ip_pb2.Constraint
Entry = ip_pb2.Entry
IP = ip_pb2.IP
Variable = ip_pb2.Variable
