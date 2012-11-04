#!/usr/bin/python3

# ***** BEGIN GPL LICENSE BLOCK *****
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ***** END GPL LICENSE BLOCK *****

# <pep8 compliant>

# Update "languages" text file used by Blender at runtime to build translations menu.

import os
import sys
import shutil

try:
    import settings
    import utils
except:
    from . import (settings, utils)

TRUNK_PO_DIR = settings.TRUNK_PO_DIR
TRUNK_MO_DIR = settings.TRUNK_MO_DIR

LANGUAGES_CATEGORIES = settings.LANGUAGES_CATEGORIES
LANGUAGES = settings.LANGUAGES
LANGUAGES_FILE = settings.LANGUAGES_FILE

OK = 0
MISSING = 1
TOOLOW = 2
FORBIDDEN = 3
FLAG_MESSAGES = {
    OK: "",
    MISSING: "No translation yet!",
    TOOLOW: "Not enough advanced to be included...",
    FORBIDDEN: "Explicitly forbidden!",
}

def find_matching_po(languages, stats, forbidden):
    """Match languages defined in LANGUAGES setting to relevant po, if possible!"""
    ret = []
    for uid, label, org_key, long_loc in languages:
        key = org_key
        if key not in stats:
            # Try to simplify the key (eg from es_ES to es).
            if '_' in org_key:
                key = org_key[0:org_key.index('_')]
            # For stuff like sr_SR@latin -> sr@latin...
            if '@' in org_key:
                key = key + org_key[org_key.index('@'):]
        if key in stats:
            if key in forbidden:
                ret.append((stats[key], uid, label, org_key, long_loc, FORBIDDEN))
            else:
                ret.append((stats[key], uid, label, org_key, long_loc, OK))
        else:
            ret.append((0.0, uid, label, org_key, long_loc, MISSING))
    return ret

def main():
    import argparse
    parser = argparse.ArgumentParser(description=""
                        "Update 'languages' text file used by Blender at runtime to build translations menu.")
    parser.add_argument('-m', '--min_translation', type=int, default=-100,
                        help="Minimum level of translation, as a percentage "
                             "(translations below this are commented out in menu).")
    parser.add_argument('langs', metavar='ISO_code', nargs='*',
                        help="Unconditionally exclude those languages from the menu.")
    args = parser.parse_args()

    ret = 0
    min_trans = args.min_translation / 100.0
    forbidden = set(args.langs)
    # 'DEFAULT' and en_US are always valid, fully-translated "languages"!
    stats = {"DEFAULT": 1.0, "en_US": 1.0}

    # Get the "done level" of each po in trunk...
    for po in os.listdir(TRUNK_PO_DIR):
        if po.endswith(".po") and not po.endswith("_raw.po"):
            lang = os.path.basename(po)[:-3]
            u1, u2, _stats = utils.parse_messages(os.path.join(TRUNK_PO_DIR, po))
            stats[lang] = _stats["trans_msg"] / _stats["tot_msg"]

    # Generate languages file used by Blender's i18n system.
    # First, match all entries in LANGUAGES to a lang in stats, if possible!
    stats = find_matching_po(LANGUAGES, stats, forbidden)
    limits = sorted(LANGUAGES_CATEGORIES, key=lambda it: it[0], reverse=True)
    idx = 0
    stats = sorted(stats, key=lambda it: it[0], reverse=True)
    langs_cats = [[] for i in range(len(limits))]
    highest_uid = 0
    for prop, uid, label, key, long_loc, flag in stats:
        if prop < limits[idx][0]:
            # Sub-sort languages by iso-codes.
            langs_cats[idx].sort(key=lambda it: it[2])
            idx += 1
        if prop < min_trans and flag == OK:
            flag = TOOLOW
        langs_cats[idx].append((uid, label, key, long_loc, flag))
        if abs(uid) > highest_uid:
            highest_uid = abs(uid)
    # Sub-sort last group of languages by iso-codes!
    langs_cats[idx].sort(key=lambda it: it[2])
    with open(os.path.join(TRUNK_MO_DIR, LANGUAGES_FILE), 'w', encoding="utf-8") as f:
        f.write("# File used by Blender to know which languages (translations) are available, \n")
        f.write("# and to generate translation menu.\n")
        f.write("#\n")
        f.write("# File format:\n")
        f.write("# ID:MENULABEL:ISOCODE:WINCODE\n")
        f.write("# ID must be unique, except for 0 value (marks categories for menu).\n")
        f.write("# Line starting with a # are comments!\n")
        f.write("#\n")
        f.write("# Automatically generated by bl_i18n_utils/update_languages_menu.py script.\n")
        f.write("# Highest ID currently in use: {}\n".format(highest_uid))
        for cat, langs_cat in zip(limits, langs_cats):
            f.write("#\n")
            # Write "category menu label"...
            if langs_cat:
                f.write("0:{}::\n".format(cat[1]))
            else:
                # Do not write the category if it has no language!
                f.write("# Void category! #0:{}:\n".format(cat[1]))
            # ...and all matching language entries!
            for uid, label, key, long_loc, flag in langs_cat:
                if flag == OK:
                    f.write("{}:{}:{}:{}\n".format(uid, label, key, long_loc))
                else:
                    # Non-existing, commented entry!
                    f.write("# {} #{}:{}:{}:{}\n".format(FLAG_MESSAGES[flag], uid, label, key, long_loc))


if __name__ == "__main__":
    print("\n\n *** Running {} *** \n".format(__file__))
    sys.exit(main())
