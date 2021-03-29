"""Send an email via Amazon SES or SMTP."""
from __future__ import print_function

from email.mime.text import MIMEText
from email.utils import formatdate
from smtplib import SMTP, SMTP_SSL

from wphase.aws import send_email_via_ses
from wphase import logger

def send_email(recipients, subject, message, from_email=None, method='ses', **kwargs):
    """Send a simple HTML email via either SES or SMTP."""
    if method == 'ses':
        send = send_email_via_ses
    elif method == 'smtp':
        send = send_email_via_smtp
    else:
        raise ValueError('Unknown email method {}'.format(method))

    if isinstance(recipients, list):
        emails = recipients
    else:
        emails = [recipients]

    if not emails:
        return

    if from_email is None:
        from_email = emails[0] # TODO is this REALLY sensible?

    # create the message
    msg = MIMEText(message, 'html')

    # add headers
    msg["From"] = from_email
    msg["To"] = ','.join(emails)
    msg["Date"] = str(formatdate(localtime=True))
    msg["Subject"] = subject

    send(recipients, msg, from_email=from_email, **kwargs)


def send_email_via_smtp(recipients, payload,
                        from_email,
                        server, port, ssl=False, tls=False,
                        user=None, password=None,
                        **kwargs):
    """Send an email (provided as a MIME payload) via SMTP."""
    if ssl:
        server = SMTP_SSL(server, port)
    else:
        server = SMTP(server, port)
        if tls:
            logger.info("SMTP: starting TLS session")
            resp, msg = server.starttls()
            if resp != 220:
                raise Exception("SMTP: could not enable TLS: %s" % msg)

    if user:
        logger.info("SMTP: attempting login")
        resp, msg = server.login(user, password)
        if resp != 235:
            raise Exception("SMTP: could not login: %s" % msg)

    refused = server.sendmail(from_email, recipients, payload.as_string())
    if refused:
        # log refused mail addresses (if any), which are reported as map of
        # form: { 'address', (error code, 'message') }
        for address, val in refused.items():
            logger.error("SMTP: %s refused by server (error %d): %s", (address, val[0], val[1]))

    resp, msg = server.quit()
    if resp != 221:
        logger.warning("SMTP: could not logout properly: %s" % msg)
