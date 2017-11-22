import obspy
from glob import glob

dirs = glob('*')
for dirname in dirs:
    filenames = glob(dirname + '/' + '*.BHZ.SAC')
    for filename in filenames:
        st = obspy.read(filename,headonly=True)
        start = st[0].stats.starttime
        common = "%d.%03d.%02d.%02d.%02d.00"%(
        start.year,start.julday,start.hour,start.minute,start.second)
        if common != filename[15:-8]:
            print(common)
            print(filename)
