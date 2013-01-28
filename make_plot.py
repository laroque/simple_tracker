#/usr/bin/python

energies = [17999 + 0.1 * i for i in range(21)]
t_final = 1e5

asim = Popen(['./aspen', str(t_final), str(energies[0])], stdout=PIPE)
asim.wait()
lines = [line.decode().strip() for line in asim.stdout.readlines()]
