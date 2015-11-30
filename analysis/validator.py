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
    parser.add_argument('trackfile', metavar='TRACKFILE', type=str, help='track file to be validated')
    parser.add_argument('-e', '--event', metavar='EVENTFILE', type=str, required=True, help='Event dat file including MCP table')
    parser.add_argument('-b', '--bintracks', action="store_true", default=True, help='treat track file as binary (default True)')
    parser.add_argument('-v', '--verbose', dest="verbose", action="store_true", default=False)
    args = parser.parse_args()

    verbose = args.verbose


    if not os.path.exists(args.trackfile):
        print("File doesn't exist: %s." % args.trackfile)
        exit(errno.ENOENT)
    if not os.path.exists(args.event):
        print("File doesn't exist: %s." % args.eventfile)
        exit(errno.ENOENT)

    event = event_model.read_datfile(args.event)
    if args.bintracks:
        tracks = event_model.read_bin_trackfile(args.trackfile, event)
    else:
        tracks = event_model.read_txt_trackfile(args.trackfile, event)
    tracks = list(tracks) # here we prefer tracks to be a list.

    weights = comp_weights(tracks, event.particles, event.hit_to_mcp)


if __name__ == "__main__":
    main()
