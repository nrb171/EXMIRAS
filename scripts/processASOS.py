import datetime
import os
import sys
from metpy.io import parse_metar_to_dataframe
from pandas import DataFrame, concat


# This script processes ASOS data files and converts them into a CSV format.
def processASOS(ASOSdir, ASOSfile, ASOSout):
    with open(ASOSdir + ASOSfile, 'r') as f:
        lines = f.readlines()

    df = DataFrame()
    for i, line in enumerate(lines):
        lineToSave = []
        timeString = line[13:27]
        # reformat the line to be parsable by Metar module
        stripped = 'METAR '+line[53:]

        dfTemp = parse_metar_to_dataframe(stripped)

        dfTemp['date_time'] = datetime.datetime.strptime(timeString, '%Y%m%d%H%M%S')

        df = concat([df,dfTemp])
        

    df.to_csv(ASOSout, index=False)


# example usage
processASOS('/h/eol/nbarron/work/asos/', '64010KOKC201506.dat', './evapModel/KOKC201506.csv')


