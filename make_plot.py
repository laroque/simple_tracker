#/usr/bin/python

#standard
from subprocess import Popen, PIPE
#3rd party
#internal
import pypeline

energies = [17999 + 0.1 * i for i in range(21)]
t_final = 5e5

asim = Popen(['./aspen', str(t_final), str(18001)], stdout=PIPE)
asim.wait()
lines = [line.decode().strip() for line in asim.stdout.readlines()]
energies = []
dphases = []
for line in lines:
    if line.split()[0] == 'time':
        energies.append(line.split()[7])
        dphases.append(line.split()[8])
datasets = sorted(zip(dphases, energies))

#the part to make the plot
plot = pypeline.usegnuplot.Gnuplot()
plot.gp("set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2")
plot.gp("set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2")
plot.gp("set style line 11 lc rgb '#808080' lt 1")
plot.gp("set border 3 back ls 11")
plot.gp("set tics nomirror")
plot.gp("set style line 12 lc rgb '#808080' lt 0 lw 1")
plot.gp("set grid back ls 12")
plot.gp("set xlabel \"dphase\"")
plot.gp("set ylabel \"ene\"")
plot.gp("set yrange [18000:18002]")
plot.gp("unset key")
plot.plot1d(datasets, ' with lines')
raw_input('waiting for you to finish looking')
