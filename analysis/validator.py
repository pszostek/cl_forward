#!/usr/bin/env python

import argparse
import errno
import os
import event_model
import numpy as np
import itertools

def comp_weights(tracks, particles, hit_to_mcp):
    w = np.zeros((len(tracks), len(particles)))
    for i, j in itertools.product(range(len(tracks)), range(len(particles))):
        trackhits = tracks[i].hits
        nhits = len(trackhits)
        particle = particles[j]
        nhits_from_p = len([h for h in trackhits if hit_to_mcp[h].count(particle) > 0])
        w[i,j] = float(nhits_from_p)/nhits
    return w

def track_association(tracks, particles, weights):
    track_to_particle = {}
    for i in range(len(tracks)):
        wtp, nwtp = np.max(weights[i,:]), np.argmax(weights[i,:])
        if wtp > 0.7:
            track_to_particle[tracks[i]] = particles[nwtp]
    return track_to_particle



def main():
    """The main function"""
    parser = argparse.ArgumentParser()
    parser.add_argument('trackfiles', metavar='TRACKFILE', nargs='+', type=str, help='txt or bin track file to be validated')
    parser.add_argument('-e', '--events', metavar='EVENTFILE', nargs='+', type=str, help='event file to be used in conjunction with txt track file')
    parser.add_argument('-b', '--bintracks', action="store_true", default=True, help='treat track file as binary (default True)')
    parser.add_argument('-v', '--verbose', dest="verbose", action="store_true", default=False)
    args = parser.parse_args()

    verbose = args.verbose

    for trackfile in args.trackfiles:
        if not os.path.exists(trackfile):
            print("Track file doesn't exist: %s." % trackfile)
            exit(errno.ENOENT)

    tracking_data = []

    if args.bintracks:
        for trackfile in args.trackfiles:
            event, tracks = event_model.read_bin_trackfile(trackfile)
            tracks = list(tracks) # here we prefer tracks to be a list.
            tracking_data.append((event, tracks))
    else:
        if len(args.trackfiles) != len(args.events):
            print("Error: Same number of event and txt track files must be provided.")
            exit(errno.ENOENT)
        for trackfile, eventfile in zip(args.trackfiles, args.events):
            if not os.path.exists(eventfile):
                print("Evnet file doesn't exist: %s." % args.trackfile)
                exit(errno.ENOENT)
            event = event_model.read_datfile(eventfile)
            tracks = event_model.read_txt_trackfile(trackfile, event)
            tracks = list(tracks) # here we prefer tracks to be a list.
            tracking_data.append((event, tracks))

    for e,ts in tracking_data:
        for t in ts:
            print t
    #weights = comp_weights(tracks, event.particles, event.hit_to_mcp)
    #track_to_particle = track_association(tracks, event.particles, weights)

# calcualte hit efficiency


if __name__ == "__main__":
    main()
