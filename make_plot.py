#/usr/bin/python2

#standard
from subprocess import Popen, PIPE
from os import getcwd
#3rd party
#internal

# global ish things
file_prefix = getcwd()
energies = [17999 + 0.1 * i for i in range(21)]
t_final = 1e4

print('foop')
filenames = []
for energy in energies:
    asim = Popen(['./aspen', str(t_final), str(18001)], stdout=PIPE)
    asim.wait()
    lines = [line.rstrip() for line in asim.stdout.readlines() if line.startswith('Finished')]
    filenames.append([file_prefix + '/' + line.split()[-1] for line in lines])
filenames = [item for sublist in filenames for item in sublist]
print('doop')

#fix from here on

g=Popen("gnuplot", stdin=PIPE, shell=False)
# do some formatting stuff here
g.stdin.write('set title "Title"\n')
g.stdin.write('set tics nomirror\n')
g.stdin.write('set xlabel \"dphase\"\n')
g.stdin.write('set ylabel \"energy\"\n')
# add the data here
file_cmd_list = []
for filename in filenames:
    (time,i,energy) = str(filename).split('.txt')[0].split('timeF')[-1].partition('energyI')
    formatstr = " using 1:2 title 'Energy " + energy + "' with lines"
    file_cmd_list.append("'" + str(filename) + "'" + formatstr)
g.stdin.write('plot ' + ', \\\n'.join(file_cmd_list) + '\n')
raw_input('waiting for you to finish looking')

#lines = [line.decode().strip() for line in asim.stdout.readlines()]
#energies = []
#dphases = []
#print('loopboop')
#for line in lines:
#    if line.split()[0] == 'time':
#        energies.append(line.split()[7])
#        dphases.append(line.split()[8])
#datasets = sorted(zip(dphases, energies))
#print('fin')
#
##the part to make the plot
#plot = usegnuplot.Gnuplot()
#plot.gp("set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2")
#plot.gp("set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2")
#plot.gp("set style line 11 lc rgb '#808080' lt 1")
#plot.gp("set border 3 back ls 11")
#plot.gp("set tics nomirror")
#plot.gp("set style line 12 lc rgb '#808080' lt 0 lw 1")
#plot.gp("set grid back ls 12")
#plot.gp("set xlabel \"dphase\"")
#plot.gp("set ylabel \"ene\"")
#plot.gp("set yrange [18000:18002]")
#plot.gp("unset key")
#plot.plot1d(datasets, ' with lines')
