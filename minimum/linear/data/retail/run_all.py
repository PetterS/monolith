import os
import glob
import subprocess


def main():
	output_file = open("output_run_all.log", "w")
	for filename in glob.glob(os.path.join(os.path.dirname(__file__),
	                                       "*.txt")):
		base_name = filename.rsplit(".", 1)[0]
		print(base_name)
		result = subprocess.run(["retail_scheduling_local", base_name],
		                        stdout=subprocess.PIPE)
		result.check_returncode()
		name, value, time = result.stdout.decode("utf-8").strip().split(" ")
		n = int(name.split("_")[0])
		output_file.write(f"{n:04d} {name} {value} {time}\n")


if __name__ == "__main__":
	main()
