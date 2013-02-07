import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import fileinput
import re
import sys

def init():
    """ Initialize the animation """


def update_animation(i):
    """ Generate one frame of animation. """
    global times, x, y, z, id_to_colour, points, ax, ttl
    t = times[i]

    to_plot_x = []
    to_plot_y = []
    to_plot_z = []
    to_plot_colour = []

    for myid in x[i].keys():
        to_plot_x.append(x[i][myid])
        to_plot_y.append(y[i][myid])
        to_plot_z.append(z[i][myid])
        to_plot_colour.append(id_to_colour[myid])

    for i in range(0, n):
        points[i].set_data(to_plot_x[i], to_plot_y[i])
        points[i].set_3d_properties(zs=to_plot_z[i])
        points[i].set_color(to_plot_colour[i])
    plt.title("Time = %f" % t)

    return points

def findmin(q):
    mins = []
    for p1 in q:
        for p2 in p1.values():
            mins.append(p2)
    return min(mins)

def findmax(q):
    maxs = []
    for p1 in q:
        for p2 in p1.values():
            maxs.append(p2)
    return max(maxs)

in_block = False
colours = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
id_to_colour = {}
times = []
x = []
y = []
z = []

coordre = re.compile("### id=(\d+), x=(-?\d+\.?\d*), y=(-?\d+\.?\d*), z=(-?\d+\.?\d*)")

for line in fileinput.input():
    if re.search("### START ENCOUNTER ###", line):
        in_block = True
    elif re.search("### END ENCOUNTER ###", line):
        in_block = False
        break
    elif re.search("### snapshot at time (\d+\.?\d*)", line):
        times.append(float(re.search("### snapshot at time (\d+\.?\d*)",
                                     line).group(1)))
        x.append({})
        y.append({})
        z.append({})
    elif coordre.search(line):
        myid = int(coordre.search(line).group(1))
        x[-1][myid] = float(coordre.search(line).group(2))
        y[-1][myid] = float(coordre.search(line).group(3))
        z[-1][myid] = float(coordre.search(line).group(4))

if len(times) == 0:
    print "Failed to find a close encounter in this data."
    sys.exit(1)

n = len(x[0].keys())

for i in range(0, n):
    id_to_colour[x[0].keys()[i]] = colours[i % len(colours)]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', 
                     xlim=(findmin(x), findmax(x)), 
                     ylim=(findmin(y), findmax(y)),
                     zlim=(findmin(z), findmax(z)))
ax.autoscale(False)

points = [None]*n
for i in range(0,n):
    points[i], = ax.plot([], [], [], 'bo')
plt.title("Time = %f" % times[0])

ani = animation.FuncAnimation(fig, update_animation, frames=len(times),
                              interval=10)

plt.show()

