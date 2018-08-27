import itertools
import time

from minimum.linear.ip import IP
from minimum.linear.solver import Solver
from minimum.linear.sum import Sum

num_groups = 8
group_size = 4

num_players = num_groups * group_size
all_players = list(range(num_players))
all_groups = list(itertools.combinations(all_players, group_size))
print("There are", num_players, "players and", len(all_groups),
      "possible groups.")


def find_groups(num_weeks):
	print("Number of weeks:", num_weeks)
	all_weeks = list(range(num_weeks))

	start_time = time.time()

	ip = IP()
	x = ip.add_boolean_grid(len(all_groups), num_weeks)

	# Symmetry breaking
	start = set()
	for i in range(num_groups):
		start.add(tuple(range(group_size * i, group_size * (i + 1))))
	second_week = tuple([group_size * i for i in range(group_size)])
	for i, g in enumerate(all_groups):
		# First week is completely fixed.
		if g in start:
			ip.add_constraint(x[i, 0] == 1)
			for w in range(1, num_weeks):
				ip.add_constraint(x[i, w] == 0)
		else:
			ip.add_constraint(x[i, 0] == 0)

		# Second week 0 and 5 plays together.
		if num_weeks > 1:
			if second_week == g:
				ip.add_constraint(x[i, 1] == 1)
			else:
				for j in second_week:
					if j in g:
						ip.add_constraint(x[i, 1] == 0)

	# Everyone plays exactly once every week.
	for w in all_weeks:
		for p in all_players:
			s = Sum()
			for i, g in enumerate(all_groups):
				if p in g:
					s += x[i, w]
			ip.add_constraint(s == 1)

	# Everyone plays together no more than once.
	for p1, p2 in itertools.combinations(all_players, 2):
		s = Sum()
		for i, g in enumerate(all_groups):
			if p1 in g and p2 in g:
				for w in all_weeks:
					s += x[i, w]
		ip.add_constraint(s <= 1)

	print("Creation time:", int(time.time() - start_time))
	start_time = time.time()

	if not Solver().solutions(ip).get():
		return False

	print("Solve time:", int(time.time() - start_time))

	for w in all_weeks:
		print("Week", w + 1)
		for i, g in enumerate(all_groups):
			if x[i, w].value() > 0.5:
				print(g, end=" ")
		print("")

	return True


w = 1
while find_groups(w):
	w += 1
