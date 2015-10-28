#!/usr/bin/env python





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
    """A simple cl_forward output file parser."""
    tracks = []
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
                tracks.append(track)
            line = tf.readline()
    return tracks

def main():
    """The main functiln"""
    trackfile1 = "results_baseline/0.out"
    tracks1 = read_trackfile(trackfile1)
    trackfile2 = "bin/x86_64/Debug/results/0_serial.out"
    tracks2 = read_trackfile(trackfile2)

    for i,t2 in enumerate(tracks2):
        eq = False
        for k,t1 in enumerate(tracks1):
            eq = eq or (t1 == t2)
        if not eq:
            print(t2)

if __name__ == "__main__":
    main()
