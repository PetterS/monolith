from IPython.core.display import HTML

from minimum.linear.variable import Variable
from minimum.linear.sum import Sum


def display_2d(x, format="%.0f", color=None, row_name_function=None):
	output = "<table>"
	output += "<tr>"
	output += "<td style='border-style:hidden;'></td>"
	for t in range(len(x[0])):
		s = str(t)
		for i in range(3 - len(s)):
			s = "&ensp;" + s
		output += "<td style='border-style:hidden;'>" + s + "</td>"
	output += "</tr>"

	for i in range(len(x)):
		output += "<tr>"
		output += "<td style='border-style:hidden;'>"
		if row_name_function is None:
			output += str(i)
		else:
			output += str(row_name_function(i))
		output + "</td>"
		for t in range(len(x[i])):
			if isinstance(x[i][t], Variable) or isinstance(x[i][t], Sum):
				v = x[i][t].value()
			else:
				v = x[i][t]

			if color is None:
				marked = abs(v) > 1e-9
			else:
				marked = color[i][t] != 0

			if marked:
				output += "<td bgcolor=\"#ccccff\" style='text-align: center;'>"
			else:
				output += "<td bgcolor=\"#eeeeee\" style='text-align: center;'>"

			if abs(v) > 1e-9:
				output += format % v
			output += "</td>"
		output += "</tr>"
	output += "</table>"
	return HTML(output)
