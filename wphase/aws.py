import os
import json
import urllib
from email.Utils import formatdate
from email.mime.text import MIMEText
import boto3

WPHASE_WEB_BASE_PATH = 'https://s3-ap-southeast-2.amazonaws.com/wphase.gagempa.net/web/'

def write_to_s3(
        output_dir,
        bucket,
        evid,
        postfix=None,
        extra_files=None,
        error_reporter=lambda x: None):
    """
    :param extra_files: Tuples of the form (<file-name>, <key>) which
        should be uploaded also.
    """

    def keygen(f):
        if postfix is None:
            return 'events/{}/wphase/{}'.format(evid, f)
        else:
            return 'events/{}/wphase/{}/{}'.format(evid, postfix, f)

    client = boto3.client('s3')
    for root, dirs, files in os.walk(output_dir):
        for f in files:
            try:
                client.upload_file(
                    os.path.join(output_dir, root, f),
                    bucket,
                    keygen(f))
            except Exception as e:
                error_reporter(str(e))

    if extra_files is not None:
        for f, k in extra_files:
            try:
                client.upload_file(f, bucket, keygen(k))
            except Exception as e:
                error_reporter(str(e))



def send_email_via_ses(
        email_address,
        bucket_name,
        event_id,
        result_id,
        eatws_env,
        call_succeeded,
        email_aws_region = 'us-west-2',
        from_email=None):

    """
    Send an email using AWS SES.

    :param str email_address: The email address to send the email to, OR a list
        of such addresses.
    :param str email_aws_region: The name of the region to use for sending the email.
        Note that SES is not available in all regions and that this will probably
        be different to the parameter *aws_region* in
        :py:func:`get_credentials`.
    """

    if isinstance(email_address, list):
        emails = email_address
    else:
        emails = [email_address]

    if not emails:
        return

    if from_email is None:
        from_email = emails[0] # TODO is this REALLY sensible?

    client = boto3.client('ses', region_name=email_aws_region)

    # all arguments to pass to the link
    all_args = {
        'success': 'true' if call_succeeded else 'false',
        'event_id': event_id,
        'result_id': result_id,
        'bucket_name': bucket_name}

    # the full url to the page
    url = WPHASE_WEB_BASE_PATH + 'index.html?' + urllib.urlencode(all_args)

    # create the message (just a link to the page to view)
    msg = MIMEText(u'<a href="{}">Go to result</a>'.format(url), 'html')

    # add headers
    msg["From"] = from_email
    msg["To"] = ','.join(emails)
    msg["Date"] = str(formatdate(localtime=True))
    msg["Subject"] = '{}W-Phase result {} ({})'.format(
        '[TEST] ' if eatws_env.lower() != 'prod' else '',
        event_id,
        result_id)

    # send the email.
    client.send_raw_email(
        Destinations = emails,
        RawMessage={'Data': msg.as_string()})
