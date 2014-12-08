#!/usr/bin/env python

"""
Simply parses the GPS data. Removes duplicates, entries with
zero velocity and accuracy less than a threshold in meters.

Extracts z(t) with t in seconds and z in meters.
"""

# threshold on accuracy
accuracy_thr = 4

file_path = './time_altitude_velocity_raw_data.dat'

from os import path
import numpy as np
import csv

unique_time_steps = []

def parse_data():

    assert path.exists(file_path)
    assert path.isfile(file_path)

    with open(file_path, 'r') as fin:
        csvr = csv.reader(fin)
        rows = [row for row in csvr]

    time, altitude, velocity = [], [], []
    duplicates = [rows[0]]

    for row in rows:
        if row in duplicates:
            continue
        else:
            duplicates += row
            t, zt, speed, accuracy = row[1], row[4], row[5], row[6]
            if t in unique_time_steps:
                continue
            else:
                unique_time_steps.append(t)
                if float(speed) <= 0. or int(accuracy) < accuracy_thr:
                    continue
                else:
                    time += [float(t)/1000.]
                    altitude += [float(zt)]
                    velocity += [float(speed)]

    # shifting data to have a zero time and zero altitude at start
    time = np.array(time, dtype=np.float) - time[0]
    altitude = np.array(altitude, dtype=np.float) - altitude[0]
    velocity = np.array(velocity, dtype=np.float)

    # removing flat part of altitude curve after 5395 seconds
    mask = (time <= 5395.)
    time = time[mask]
    altitude = altitude[mask]
    velocity = velocity[mask]

    # getting the first part of the graph (ascent)
    ascent_time = (0 <= time) * (time <= 2880)
    ascent = altitude[ascent_time]

    # getting the second part of the graph (descent)
    descent_time = (4160 <= time) * (time <= 5250)
    descent = altitude[descent_time]

    return ascent, descent, time[ascent_time], time[descent_time]


def main():
    # plot resulting altitude with time
    import matplotlib.pyplot as plt

    ascent, descent, ascent_time, descent_time = parse_data()

    plt.plot(ascent_time, ascent, 'r-', label='ascent')
    plt.plot(descent_time, descent, 'b-', label='descent')
    plt.xlabel('time (s)')
    plt.ylabel('altitude (m)')
    plt.show()


if __name__ == '__main__':
    main()
