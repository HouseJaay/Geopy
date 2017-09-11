import datetime
import os
from email.mime.text import MIMEText
import smtplib
from time import sleep
import random
import sys
import haotool as ht

sys.path.insert(0,'/home/hao_shijie/work/two_station')

import cal_station as cs

def def_mail(row_pair,evt):
    with open('config','r') as f:
        header = f.read()
    label = row_pair['station1'] +'_'+ row_pair['station2']
    header = header + '.LABEL ' + label + '\n.QUALITY B\n.END\n'
    def def_window(evtime):
        begin = evtime
        end = evtime + 5000
        data = " %d %02d %02d %02d %02d 00.0 %d %02d %02d %02d %02d 00.0 1 BHZ\n" \
        % (begin.year,begin.month,begin.day,begin.hour,begin.minute,\
        end.year,end.month,end.day,end.hour,end.minute)
        return data
    
    body = ''
    sta1 = row_pair['station1'] +' '+ row_pair['net1']
    sta2 = row_pair['station2'] +' '+ row_pair['net2']
    for i in evt.index:
        evtime = evt.loc[i]['time']
        req = def_window(evtime)
        body += sta1 + req
        body += sta2 + req
    return header+body

def send_mail(text):
    msg = MIMEText(text,'plain','utf-8')
    from_addr = "seishao@126.com"
    #password = "xiewanlehaha"
    password = "fuckyou250"
    #to_addr = "xlsthsj@126.com"
    to_addr = "breq_fast@iris.washington.edu"
    smtp_server = "smtp.126.com"
    msg['From'] = from_addr
    msg['To'] = to_addr
    msg['Subject'] = "Thesis task"
    server = smtplib.SMTP(smtp_server, 25)
    server.set_debuglevel(1)
    server.login(from_addr, password)
    server.sendmail(from_addr, [to_addr], msg.as_string())
    server.quit()

def do_send(pair,evt_full):
    for i in pair.index:
        row_pair = pair.loc[i]
        evt = cs.get_event(row_pair,evt_full,dep_max,dist_min,dist_max,mag_min)
        print(len(evt))
        text = def_mail(row_pair,evt)
        send_mail(text)
        print('send mail for data %s %s' %(row_pair['station1'],row_pair['station2']))
        sleep(random.randrange(30,90))

dep_max = 30
dist_min = 2000
dist_max = 9000
mag_min = 5.8

if __name__=="__main__":
    evt_full = ht.read_event('evt_2004_2016')
    sta = ht.read_station('station_us')
    pair = ht.mk_sta_pairs(sta)
    pair.loc[:,'event'] = pair.apply(cs.do_check,axis='columns',
    args=(evt_full,dep_max,dist_min,dist_max,mag_min))
    pair_temp = pair[pair['event']>5]
    #do_send(pair_temp,evt_full)
