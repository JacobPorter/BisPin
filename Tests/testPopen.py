import subprocess
myp = subprocess.Popen(["cat", "/proc/cpuinfo"], stdout = subprocess.PIPE)
myp2 = subprocess.Popen(["wc", "-l"], stdin = myp.stdout)
myp2.communicate()
#for line in myp.stdout:
#	myp2.communicate(line)

