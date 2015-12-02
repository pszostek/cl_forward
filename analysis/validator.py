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

def hit_purity(tracks, particles, weights):
    track_to_particle = {t:(0.0, None) for t in tracks}
    for i in range(len(tracks)):
        wtp, nwtp = np.max(weights[i,:]), np.argmax(weights[i,:])
        if wtp > 0.7:
            track_to_particle[tracks[i]] = (wtp, particles[nwtp])
        else:
            track_to_particle[tracks[i]] = (wtp, None)
    return track_to_particle

def hit_efficinecy(track_to_particle, hit_to_mcp, mcp_to_hits):
    hit_eff = {}
    for track, (hp, particle) in track_to_particle.iteritems():
        if particle is None:
            continue # no particle assoicated
        hits_p_on_t = sum([hit_to_mcp[h].count(particle) for h in track.hits])
        hit_eff[(track, particle)] = float(hits_p_on_t)/len(mcp_to_hits[particle])
    return hit_eff

def track_recoeff(track_to_particle, particles):
    nreconstructible = len(particles)
    nreconstructed = len(set([pp[1] for _,pp in track_to_particle.iteritems() if pp[1] is not None]))
    return float(nreconstructed)/nreconstructible

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
    if verbose:
        print("Reading data: ")
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

    if verbose:
        print(" done.")
    avg_recoeff = 0.0
    for i, (event, tracks) in enumerate(tracking_data):
        weights = comp_weights(tracks, event.particles, event.hit_to_mcp)
        track_to_particle = hit_purity(tracks, event.particles, weights)
        hit_eff = hit_efficinecy(track_to_particle, event.hit_to_mcp, event.mcp_to_hits)
        #for k,v in hit_eff.iteritems():
        #    print k[0].tid,k[1].pkey,v
        recoeff = track_recoeff(track_to_particle, event.particles)
        print recoeff
        avg_recoeff += recoeff
    print("Average reconstruction efficiency: ",avg_recoeff/len(tracking_data))


if __name__ == "__main__":
    main()
