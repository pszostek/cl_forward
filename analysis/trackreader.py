#!/usr/bin/env python

import argparse
import errno
import os


class Hit(object):
    """Hit instances store all relevant information of hits from the track"""

    def __init__(self, hitline):
        """Reads a hit line from the cl_forward output file and creates a Hit object"""
        fields = hitline.strip().split()
        self.hitID = int(fields[0])
        self.hitNum = int(fields[1][1:-1])
        self.x = float(fields[5][:-1])
        self.y = float(fields[7][:-1])
        self.z = float(fields[9])

    def __str__(self):
        return "%d (%d): x = %g, y = %g, z = %g"%(self.hitID, self.hitNum, self.x, self.y, self.z)

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return self.hitID


class Track(object):
    """Track instances hold all information needed to reconstruct a track, that
    includes also all the hits that make up this track"""

    def __init__(self, tid, lines):
        """Takes a section of cl_forward output file that contains all the track hits"""
        self.tid = tid
        self.hits = []
        for line in lines:
            hit = Hit(line)
            self.hits.append(hit)

    def get_coords(self):
        """Return three lists for x, y and z coordinates of hits.
        This is especially useful for using with matplotlib 3D plotting."""
        xs = []
        ys = []
        zs = []
        for hit in self.hits:
            xs.append(hit.x)
            ys.append(hit.y)
            zs.append(hit.z)
        return (xs, ys, zs)

    def __str__(self):
        s = 'Track (#%d) length: %d'%(self.tid,len(self.hits))
        for hit in self.hits:
            s += '\n' + str(hit)
        return s

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if len(self.hits) != len(other.hits):
                return False
            eq = True
            for i, h in enumerate(self.hits):
                eq = eq and (h == other.hits[i])
            return eq
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        hsh = 0
        for i,h in enumerate(self.hits):
            hsh += (i+1)*h.hitID
        return hsh

def read_trackfile(filename):
    """A simple cl_forward output file parser.

    returns a set of tracks in the file
    """
    tracks = set()
    with open(filename) as tf:
        line = tf.readline()
        while line:
            if line.startswith('Track'):
                tid = int(line.strip().split()[1][1:-1])
                tlen = int(line.strip().split()[-1])
                tracklines = []
                for _ in xrange(tlen):
                    tracklines.append(tf.readline())
                track = Track(tid,tracklines)
                tracks.add(track)
            line = tf.readline()
    return tracks


def main():
    """The main function"""
    parser = argparse.ArgumentParser()
    parser.add_argument('files', metavar='FILE', type=str, nargs='+', help='files to be compared')
    parser.add_argument('--print-difference', dest='print_difference', action="store_true", default=False, help='Print tracks present in the second file,'
                        'but not in the first one')
    parser.add_argument('--test-equal', dest="test_equal", action="store_true", default=False, help='Test two files for equality')
    parser.add_argument('--verbose', dest="verbose", action="store_true", default=False)
    parser.add_argument('--tracks', dest="tracks", action="store_true", default=False, help='Print number of tracks in the file')
    args = parser.parse_args()

    verbose = args.verbose

    for file_ in args.files:
        if not os.path.exists(file_):
            print("File doesn't exist: %s." % file_)
            exit(errno.ENOENT)

    tracks1 = read_trackfile(args.files[0])
    try:
        tracks2 = read_trackfile(args.files[1])
    except IndexError:  # it will be raised if there is only one file
        tracks2 = None

    if verbose:
        print("Number of tracks in %s: %d" % (args.files[0], len(tracks1)))
        print("Number of tracks in %s: %d" % (args.files[1], len(tracks2)))

    if args.test_equal:
        assert(tracks2 is not None)
        equal = True
        if len(tracks1) != len(tracks2):
            print("Number of tracks in two files doesn't match.")
            equal = False
        t1_minus_t2 = tracks1-tracks2
        t2_minus_t1 = tracks2-tracks1

        if len(t1_minus_t2):
            print("There are %d tracks in %s, but not in %s" % (len(t1_minus_t2), args.files[0], args.files[1]))
            if verbose:
                for track in t1_minus_t2:
                    print(track)
            equal = False

        if len(t2_minus_t1):
            print("There are %d tracks in %s, but not in %s" % (len(t2_minus_t1), args.files[1], args.files[0]))
            if verbose:
                for track in t2_minus_t1:
                    print(track)
            equal = False

        if equal is True:
            print("The two files match.")
            exit(0)
        else:
            exit(1)

    elif args.tracks:
        assert(tracks1 is not None)
        print(len(tracks1))
        exit(0)

    elif args.print_difference:

        for i,t2 in enumerate(tracks2):
            eq = False
            for k,t1 in enumerate(tracks1):
                eq = eq or (t1 == t2)
            if not eq:
                print(t2)
    else:
        print("Ooops.. How about saying what should I do?")
        exit(errno.EAGAIN)

if __name__ == "__main__":
    main()
