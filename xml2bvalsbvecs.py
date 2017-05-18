#!/usr/bin/env python3

# Converts ADNI XML files containing DTI bvalues and gradients to FSL bvals and bvecs files
# Usage:
# ./xml2bvalsbvecs bvals_file bvecs_file file1.xml [file2.xml [...]]

import re
import sys


def xmls_to_bvals_bvecs(bvals_file, bvecs_file, filenames):

    bvals = []
    xgrads = []
    ygrads = []
    zgrads = []

    for filename in filenames:
        # It would be better to parse it as XML, but the XML files are not valid
        xml = open(filename).read().replace("\n", "")

        bvals.append( float(re.match('.*bvalue=\"(.*?)\".*', xml).groups()[0]))
        xgrads.append(float(re.match('.*xgradient=\"(.*?)\".*', xml).groups()[0]))
        ygrads.append(float(re.match('.*ygradient=\"(.*?)\".*', xml).groups()[0]))
        zgrads.append(float(re.match('.*zgradient=\"(.*?)\".*', xml).groups()[0]))

    fbvals = open(bvals_file, "w")
    fbvals.write(" ".join(["%f" % val for val in bvals]) + "\n")
    fbvals.close()

    fbvecs = open(bvecs_file, "w")
    fbvecs.write(" ".join(["%f" % val for val in xgrads]) + "\n")
    fbvecs.write(" ".join(["%f" % val for val in ygrads]) + "\n")
    fbvecs.write(" ".join(["%f" % val for val in zgrads]) + "\n")
    fbvecs.close()



if __name__ == "__main__":
    xmls_to_bvals_bvecs(sys.argv[1], sys.argv[2], sys.argv[3:])