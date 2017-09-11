import datetime
import os

sta_file="sta.lst" #station list
start_time=datetime.date(2006,2,1) #time range
end_time=datetime.date(2007,11,1)
time_delta=datetime.timedelta(days=1) #file length

with open(sta_file,'r') as f:
    stas=f.read().splitlines()
for sta in stas:
    with open('config','r') as f:
        text=f.read()
    label=sta.split(' ')[1]+'_'+sta.split(' ')[0]
    text=text+'.LABEL '+label+'\n.QUALITY B\n.END\n'
    day=start_time
    while(day<end_time):
        day2=day+time_delta
        request=sta+" %d %02d %02d 00 00 00.0 %d %02d %02d 00 00 00.0 1 BHZ\n" % (day.year,day.month,day.day,day2.year,day2.month,day2.day)
        text+=request
        day=day2
    with open('mail','w') as f:
        f.write(text)
    #os.system('mailx breq_fast@iris.washington.edu < mail')
    #os.system('rm mail')
