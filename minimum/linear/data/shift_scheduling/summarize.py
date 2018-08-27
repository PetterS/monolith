#!/usr/bin/python3

import sqlite3
import glob

for file in sorted(glob.glob("*.db")):
	db = sqlite3.connect(file)
	res = db.execute(
	    "SELECT objective FROM solution ORDER BY objective ASC LIMIT 1;"
	).fetchone()
	best = res[0] if res else None
	try:
		refined = db.execute(
		    "SELECT objective FROM refined ORDER BY objective ASC LIMIT 1;"
		).fetchone()[0]
	except sqlite3.OperationalError:
		refined = None
	print(file, best, refined)
