#!/usr/bin/env python3

import os
import requests
import sys


ldraw_path = os.path.abspath(os.path.dirname(__file__))


def download_one_part(part):
    """Downloads one part or subpart."""

    # Try to read the parts file
    parts_url = "https://www.ldraw.org/library/unofficial/parts/" + part
    req = requests.get(parts_url, allow_redirects=True)
    is_part = req.content[0] == ord("0")

    # Part not found, try p file
    if not is_part:
        p_url = "https://www.ldraw.org/library/unofficial/p/" + part
        req = requests.get(p_url, allow_redirects=True)
        if req.content[0] != ord("0"):
            sys.exit("Could not find part: " + part)

    # Save the file
    local_path = os.path.join(ldraw_path, "parts" if is_part else "p", part)
    print(local_path)
    with open(local_path, "wb") as local_file:
        local_file.write(req.content)

    # Display first line, which is a description of the part
    print("Downloaded:", req.content.decode().splitlines()[0])


def get_dependencies(part):
    """Gets list of dependencies of a previously downloaded part."""
    # Try to find the part in /parts or /p
    try:
        part_contents = open(os.path.join(ldraw_path, "parts", part)).read()
    except FileNotFoundError:
        part_contents = open(os.path.join(ldraw_path, "p", part)).read()

    # Go through each line, and yield the part names, which are dependencies
    for line in part_contents.splitlines():
        *_, part_name = line.split(' ')

        # Only return part names, and reverse backslash to forward slash
        if '.dat' in part_name:
            part_name = part_name.replace("\\", "/")
            yield part_name


def download_full_part(part):
    """Downloads a part and its dependencies recursively"""

    # top level part for this iteration
    download_one_part(part)

    # Find out what it depends on, then download those too
    for dep in get_dependencies(part):

        parts_path = os.path.join(ldraw_path, "parts", dep)
        p_path = os.path.join(ldraw_path, "p", dep)

        # Only download it if it doesn't already exist
        if not os.path.isfile(parts_path) and not os.path.isfile(p_path):
            download_full_part(dep)


# Command line interface
if __name__ == "__main__":

    # Check arguments
    if len(sys.argv) != 2:
        print("Example usage:\n\n    ./download_part.py 69730c01.dat")

    # Download main part
    download_full_part(sys.argv[1])
