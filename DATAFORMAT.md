Binary file data format
=======================

We define three types of files:

-	An *event file* contains all necessary information for a velopix track reconstruction code. That is, information on sensor location and hits.
-	*Extended event files* append *event files* by optional information from the Monte-Carlo simulations, notably the MCParticles and particle-to-hit tables.
-	*Track files* are directly derived from *event files*, appending event and MC truth information with the reconstructed tracks.

We can describe the two (three) files in one table since *extended event files* and *track files* are based on simple *event files*. This design allows us to be compatible with the original data format as defined for the coprocessor manager and cl_forward.

Format description
------------------

File organization:

```

|------------------------|
|       File header      |
|------------------------|
|      Event segment     |
|------------------------|
|   MC Particle segment  |
|------------------------|
|     Track segment      |
|------------------------|



```

File header:

| name        | type     | size        |
|-------------|----------|-------------|
| funcNameLen | uint32_t | 1           |
| funcName    | char     | funcNameLen |
| dataSize    | uint32_t | 1           |

Event Segment:

| name             | type     | size         |
|------------------|----------|--------------|
| no_sensors       | uint32_t | 1            |
| no_hits          | uint32_t | 1            |
| sensor_Zs        | uint32_t | h_no_sensors |
| sensor_hitStarts | uint32_t | h_no_sensors |
| sensor_hitNums   | uint32_t | h_no_sensors |
| hit_IDs          | uint32_t | h_no_hits    |
| hit_Xs           | float    | h_no_hits    |
| hit_Ys           | float    | h_no_hits    |
| hit_Zs           | float    | h_no_hits    |

MC Particle segment:

The segment begins with one uint denoting the number of MC particles following

| name   | type     | size |
|--------|----------|------|
| no_mcp | uint32_t | 1    |

This is followed by `no_mcp` structures. One for each MC particle:

| name    | type     | size    |
|---------|----------|---------|
| mcp_key | uint32_t | 1       |
| mcp_id  | uint32_t | 1       |
| mcp_p   | float    | 1       |
| mcp_pt  | float    | 1       |
| mcp_eta | float    | 1       |
| mcp_phi | float    | 1       |
| no_hits | uint32_t | 1       |
| hitIDs  | uint32_t | no_hits |

Track segment:

Similar to the MC particle segment we begin with an integer denoting the number of tracks, followed by track structures:

| name      | type     | size |
|-----------|----------|------|
| no_tracks | uint32_t | 1    |

| name    | type     | size    |
|---------|----------|---------|
| no_hits | uint32_t | 1       |
| hitIDs  | uint32_t | no_hits |
