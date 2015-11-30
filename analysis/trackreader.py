#!/usr/bin/env python

import argparse
import errno
import os
import event_model




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

    tracks1 = event_model.read_txt_trackfile(args.files[0])
    try:
        tracks2 = event_model.read_txt_trackfile(args.files[1])
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
