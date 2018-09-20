#!/usr/bin/env python3

#*******************************************************************************
#* reducer.py
#*
#* Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
#*
#* All rights reserved. Published under the BSD-2 license in the LICENSE file.
#******************************************************************************/

import argparse
import subprocess
import copy
import random

parser = argparse.ArgumentParser()
parser.add_argument('filename')
parser.add_argument('cmd', nargs=argparse.REMAINDER)
args = parser.parse_args()

RF = "REDUCEFILE"
LRF = "LAST_REDUCEFILE"

def read_binary(path):
    with open(path, 'rb') as file:
        inp = file.read()
    return inp

def write_binary(path, data):
    with open(path, 'wb') as file:
        file.write(data)
        file.flush()
    test = read_binary(path)
    assert test == data

def run(cmd, current):
    print(cmd)
    print("Input size: {}".format(len(current)))

    write_binary(RF, current)

    print("------------------------------------")
    res = subprocess.run(cmd)
    print("------------------------------------")


    if res.returncode != 0:
        return True
    else:
        return False

inp = read_binary(args.filename)
write_binary(RF, inp)

cmd = list(args.cmd)

prev = copy.deepcopy(inp)
current = copy.deepcopy(inp)
while True:
    print("====================================")
    if run(cmd, current):
        # reduce further
        prev = copy.deepcopy(current)
        write_binary(LRF, prev)

        ri = random.randrange(0, len(current))
        rl = random.randrange(0, 100) + 1

        left = current[:ri]
        right = current[(ri+rl):]
        current = left + right
    else:
        print("Did not fail, try something different")
        if prev == current:
            print("Got stuck")
            exit(1)
        current = copy.deepcopy(prev)


