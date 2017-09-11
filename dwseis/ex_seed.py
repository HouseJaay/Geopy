from glob import glob
import os

seed_list = glob('*.seed')
for seed_file in seed_list:
    out_dir = "../" + seed_file.split('.')[0]
    command = "rdseed -R -d -z 1 -q %s -f %s" %(out_dir,seed_file)
    os.mkdir(out_dir)
    os.system(command)
