import datetime
import os
import re
import sys

from aliyunsdkcore.request import CommonRequest

from aliyunsdkcore.client import AcsClient

# edit here
remote_sms_client = AcsClient("AccessKeyID", "AccessKeySecret", 'cn-hangzhou')

def sms_remote(phone,name,job, time,total,fin,rest ):
    job=job[:20]
    time=time[:20]


    request = CommonRequest()

    request.set_accept_format('json')
    request.set_domain('dysmsapi.aliyuncs.com')
    request.set_method('POST')
    request.set_protocol_type('https')  # https | http
    request.set_version('2017-05-25')
    request.set_action_name('SendSms')
    request.add_query_param('RegionId', "cn-hangzhou")
    request.add_query_param('PhoneNumbers', phone)
    # edit here
    request.add_query_param('SignName', "SignName")
    request.add_query_param('TemplateCode', "TemplateCode")

    request.add_query_param('TemplateParam',
                                "{\"user\":\"" + name + "\",\"jobname\":\"" + job +"\",\"time\":\"" + time + "\",\"total_job_num\":\"" + total + "\",\"fin\":\"" + fin +"\",\"remain_job_num\":\"" + rest + "\"}")
    response = remote_sms_client.do_action_with_exception(request)
    print(str(response, encoding='utf-8'))


phones={
        "user1":"phone1",
        "user2":"phone2",
        "user3":"phone3",
        }

    


if __name__ == '__main__':
    who=sys.argv[1]
    job=sys.argv[2]
    total=sys.argv[3]
    fin=sys.argv[4]
    rest=sys.argv[5]
    time = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    phone=phones[who]
    sms_remote(phone,who,job,time,total,fin,rest)

