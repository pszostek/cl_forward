#!/usr/bin/env python
from struct import unpack

class Event(object):
    """A SOA datastructure for events"""

    def __init__(self, sensor_Zs, sensor_hitStarts, sensor_hitNums,
        hit_IDs, hit_Xs, hit_Ys, hit_Zs, hits, mcp_to_hits=None):
        self.sensor_Zs = sensor_Zs
        self.sensor_hitStarts = sensor_hitStarts
        self.sensor_hitNums = sensor_hitNums
        self.hit_IDs = hit_IDs
        self.hit_Xs, self.hit_Ys, self.hit_Zs = hit_Xs, hit_Ys, hit_Zs
        self.hits = hits

        self.hits_by_id = {hit.hitID:hit for hit in self.hits}

        self.mcp_to_hits = mcp_to_hits
        self.hit_to_mcp = None
        self.particles = None
        if self.mcp_to_hits is not None:
            self.particles = list(self.mcp_to_hits.keys())
            self.hit_to_mcp = {h:[] for h in self.hits}
            for mcp,mhits in mcp_to_hits.iteritems():
                for hit in mhits:
                    self.hit_to_mcp[hit].append(mcp)

    def get_hit(self, hitID):
        return self.hits_by_id[hitID]

class MCParticle(object):
    """Store information about a Monte-Carlo simulation particle"""

    def __init__(self, pkey, pid, velohits):
        """Construct a new particle from

        its numeric key (arbitrary integer used in input file)
        its pID (integer providing information on the particle type)
        its assocated velopixel hits"""
        self.pkey = pkey
        self.pid = pid
        self.velohits = velohits

    def __str__(self):
        return "MCParticle %d: pid = %d"%(self.pkey, self.pid)

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__)


class Hit(object):
    """Hit instances store all relevant information of hits from the track"""

    def __init__(self, hitID, x=0.0, y=0.0, z=0.0, module=0):
        """Construct a hit from its coordinates and ID"""
        self.hitID = hitID
        self.module = module
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return "%d M#(%d): x = %g, y = %g, z = %g"%(self.hitID, self.module, self.x, self.y, self.z)

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

    def __init__(self, tid, trackhits):
        """Constructs a new track from a list of hits"""
        self.tid = tid
        self.hits = trackhits

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
        s = 'Track (#%d) length: %d'%(self.tid, len(self.hits))
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

def read_datfile(filename):
    with open(filename, 'rb') as f:
        bindata = f.read()
    pos = 0
    # read file header
    funcNameLen, = unpack('i', bindata[pos:pos+4])
    pos += 4
    funcName = unpack('c'*funcNameLen, bindata[pos:pos+funcNameLen])
    funcName = ''.join(funcName)
    pos += funcNameLen
    dataLen, = unpack('i', bindata[pos:pos+4])
    pos += 4
    header_sz = pos
    no_sensors, no_hits = unpack('ii', bindata[pos:pos+8])
    pos += 8
    # read sensor data
    sensor_Zs = unpack('i'*no_sensors, bindata[pos:pos+no_sensors*4])
    pos += no_sensors*4
    sensor_hitStarts = unpack('i'*no_sensors, bindata[pos:pos+no_sensors*4])
    pos += no_sensors*4
    sensor_hitNums = unpack('i'*no_sensors, bindata[pos:pos+no_sensors*4])
    pos += no_sensors*4
    # read hit data
    hit_IDs = unpack('i'*no_hits, bindata[pos:pos+no_hits*4])
    pos += no_hits*4
    hit_Xs = unpack('f'*no_hits, bindata[pos:pos+no_hits*4])
    pos += no_hits*4
    hit_Ys = unpack('f'*no_hits, bindata[pos:pos+no_hits*4])
    pos += no_hits*4
    hit_Zs = unpack('f'*no_hits, bindata[pos:pos+no_hits*4])
    pos += no_hits*4
    hits = []
    for hid, x, y, z in zip(hit_IDs, hit_Xs, hit_Ys, hit_Zs):
        hits.append(Hit(hid, x, y, z))
    mcp_to_hits = {}
    if (pos < dataLen + header_sz):
        # it looks like we have some MCP data in this file, proceed
        no_mcp, = unpack('i', bindata[pos:pos+4])
        pos += 4
        for _ in range(no_mcp):
            mcpkey, mcpid, nh = unpack('iii', bindata[pos:pos+12])
            pos += 12
            mcp_hitIDs = unpack('i'*nh, bindata[pos:pos+nh*4])
            pos += nh*4
            hdict = {h.hitID:h for h in hits}
            trackhits = [hdict[hid] for hid in mcp_hitIDs]
            mcp_to_hits[MCParticle(mcpkey, mcpid, trackhits)] = trackhits
            #mcp_to_hits[mcpid] = mcp_hitIDs
    if pos < dataLen + header_sz:
        # if there is even more data, return access to the caller
        # he'll know what to do with it
        return Event(sensor_Zs, sensor_hitStarts, sensor_hitNums,
                hit_IDs, hit_Xs, hit_Ys, hit_Zs, hits, mcp_to_hits),(bindata,pos)
    else:
        return Event(sensor_Zs, sensor_hitStarts, sensor_hitNums,
                hit_IDs, hit_Xs, hit_Ys, hit_Zs, hits, mcp_to_hits)

def read_bin_trackfile(filename):
    event, (bindata, pos) = read_datfile(filename)
    nTracks, = unpack('i',bindata[pos:pos+4])
    pos += 4
    tracks = set()
    for itrack in range(nTracks):
        nhits, = unpack('i', bindata[pos:pos+4])
        pos += 4
        track_hitIDs = unpack('i'*nhits, bindata[pos:pos+nhits*4])
        if event is None:
            trackhits = [Hit(hid) for hid in track_hitIDs]
        else:
            trackhits = [event.hits_by_id[hid] for hid in track_hitIDs]
        pos += nhits*4
        tracks.add(Track(itrack, trackhits))
    return event, tracks


def read_txt_trackfile(filename, event=None):
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
                trackhits = []
                for _ in xrange(tlen):
                    hitline = tf.readline()
                    flds = hitline.strip().split()
                    hitID, module = int(flds[0]),int(flds[3][:-1])
                    #self.hitNum = int(flds[1][1:-1])
                    if event is None:
                        x, y, z = float(flds[5][:-1]), float(flds[7][:-1]), float(flds[9])
                        hit = Hit(hitID, x, y, z, module)
                    else:
                        hit = event.get_hit(hitID)
                    trackhits.append(hit)
                track = Track(tid,trackhits)
                tracks.add(track)
            line = tf.readline()
    return tracks
