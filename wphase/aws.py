from future import standard_library
standard_library.install_aliases()
from builtins import str

import json
import os
import urllib.error
import urllib.parse
import urllib.request
from email.mime.text import MIMEText
from email.Utils import formatdate

import boto3


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
        subject,
        message,
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

    # create the message
    msg = MIMEText(message, 'html')

    # add headers
    msg["From"] = from_email
    msg["To"] = ','.join(emails)
    msg["Date"] = str(formatdate(localtime=True))
    msg["Subject"] = subject

    # send the email.
    client.send_raw_email(
        Destinations = emails,
        RawMessage={'Data': msg.as_string()})
