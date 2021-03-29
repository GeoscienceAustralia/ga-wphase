#!/usr/bin/env python
"""Runs W-Phase, sends results to seiscomp3 messaging and copies files to S3."""
from __future__ import absolute_import, print_function

import logging
import os
import sys
import json

from datetime import datetime, timedelta

from obspy.core import UTCDateTime

from seiscomp3 import DataModel as DM, Logging, IO, Core
from seiscomp3.Client import Application

from wphase import runwphase, settings
from wphase.aws import write_to_s3
from wphase.email import send_email
from wphase.result import WPhaseParser, FMItem
from wphase.seiscomp import createAndSendObjects, writeSCML


# this is used in determining the name of the log file that this application
# logs to. I think this must be the name of this file, but cannot be sure as the
# logging is Gempa's black magic.
APPLICATION_NAME = 'wphase'

LOG_FILE_NAME = '{}/.seiscomp3/log/{}.log'.format(
    os.environ['HOME'], APPLICATION_NAME)

class LogHandler(logging.Handler):
    """Handler that forwards a python logger to seiscomp logs"""
    def emit(self, record):
        msg = self.format(record)
        if record.levelname == 'DEBUG':
            Logging.debug(msg.encode())
        elif record.levelname == 'INFO':
            Logging.info(msg.encode())
        elif record.levelname == 'WARNING':
            Logging.warning(msg.encode())
        else:
            Logging.error(msg.encode())

# Send W-Phase logs and python warnings to seiscomp logs
import wphase.logger
logger = wphase.logger.getLogger()
logger.setLevel(logging.DEBUG)
logger.addHandler(LogHandler())
wlogger = logging.getLogger('py.warnings')
wlogger.setLevel(logging.WARNING)
wlogger.addHandler(LogHandler())
logging.captureWarnings(True)

class WPhase(Application):
    def __init__(self, argc, argv):
        # Default to loglevel 1
        Application.__init__(self, argc, argv)

        # default location to write outputs to
        self.output = '/tmp/wphase-output'
        self.filename = None
        self.mag_type = None
        self.mag_value = None
        self.server = 'IRIS'
        self.networks = 'ALL'
        self.region = 'not specified'
        self.evid = None
        self.resultid = None
        self.notificationemail = None
        self.fromemail = None
        self.email_aws_region = None
        self.email_method = 'ses'

        self.smtp_server = None
        self.smtp_port = 25
        self.smtp_ssl = False
        self.smtp_tls = False
        self.smtp_user = None
        self.smtp_password = None

        self.write_s3 = False
        self.bucket_name = None
        self.agency = 'GA'
        self.make_maps = True

        # enable messaging support
        self.setMessagingEnabled(True)

        # disable database access
        self.setDatabaseEnabled(False, False)

        # default spread username
        self.setMessagingUsername("gawphase")

        # send all objects to the focal mechanism group
        self.setPrimaryMessagingGroup("FOCMECH")



    def createCommandLineDescription(self):
        """
        Override this method (including a super call) in a subclass to add
        custom options.
        """

        self.commandline().addGroup("Output")
        self.commandline().addStringOption(
            "Output",
            "outputs,o",
            "Directory to write output to. Defaults to /tmp/wphase-output.")

        self.commandline().addGroup("Input")
        self.commandline().addStringOption(
            "Input",
            "lat",
            "Latitude of the event.")
        self.commandline().addStringOption(
            "Input",
            "lon",
            "Longitude of the event.")
        self.commandline().addStringOption(
            "Input",
            "depth",
            "Depth of the event.")
        self.commandline().addStringOption(
            "Input",
            "sourcezone",
            "Source zone of the event.")
        self.commandline().addStringOption(
            "Input",
            "time",
            "Time of the event.")
        self.commandline().addStringOption(
            "Input",
            "magtype",
            "The type of magnitude of the triggering event.")
        self.commandline().addStringOption(
            "Input",
            "magvalue",
            "The magnitude of the triggering event.")
        self.commandline().addStringOption(
            "Input",
            "server",
            "The FDSN server to use. This can be any server supported by " + \
            "the obspy FDSN client, or an arbitrar URL.")
        self.commandline().addStringOption(
            "Input",
            "networks",
            "A comma separated list of networks to use.")
        self.commandline().addStringOption(
            "Input",
            "region",
            "The name of the region in which the event occured.")
        self.commandline().addStringOption(
            "Input",
            "evid",
            "First part of the key under which to write the events to S3.")
        self.commandline().addStringOption(
            "Input",
            "resultid",
            "The second part of the key under which to write the events to S3.")
        self.commandline().addStringOption(
            "Input",
            "notificationemail",
            "Email address to send notification to.")
        self.commandline().addStringOption(
            "Input",
            "fromemail",
            "Email address to send notification from.")
        self.commandline().addStringOption(
            "Input",
            "emailmethod",
            "Method to use to send notification email (either ses or smtp)")
        self.commandline().addStringOption(
            "Input",
            "emailawsregion",
            "If using SES, AWS Region to send email notifications from.")
        self.commandline().addStringOption(
            "Input",
            "smtp-server",
            "SMTP server to use to send notification email.")
        self.commandline().addStringOption(
            "Input",
            "smtp-port",
            "SMTP port.")
        self.commandline().addStringOption(
            "Input",
            "smtp-user",
            "SMTP username.")
        self.commandline().addStringOption(
            "Input",
            "smtp-password",
            "SMTP username.")
        self.commandline().addOption(
            "Input",
            "smtp-ssl",
            "Enable SSL for SMTP.")
        self.commandline().addOption(
            "Input",
            "smtp-tls",
            "Enable TLS for SMTP.")
        self.commandline().addStringOption(
            "Input",
            "writeS3",
            "Should results be written to S3?")
        self.commandline().addStringOption(
            "Input",
            "bucketname",
            "The name of the S3 bucket to write to.")
        self.commandline().addStringOption(
            "Input",
            "agency",
            "Agency code for SC3 output")
        self.commandline().addOption(
            "Input",
            "make-maps",
            "Whether to create map plots or not.")
        self.commandline().addOption(
            "Input",
            "overwrite",
            "Whether to overwrite existing outputs or not.")

    def processArg(self, name, to=None, default=None, conv=str):
        if to is None:
            to = name

        try:
            val = self.commandline().optionString(name)
        except Exception:
            if default is not None:
                setattr(self, to, default)
        else:
            setattr(self, to, conv(val))


    def validateParameters(self):
        """
        SeisComP3 specific method.
        """

        if not Application.validateParameters(self):
            return False

        getter = self.processArg

        try:
            # If there is an unrecognized option it must be a JSON file
            # of wphase outputs. In this case, the file is parsed and pushed
            # to the messaging system and written to disk.
            self.filename = self.commandline().unrecognizedOptions()[0]
        except Exception:
            # Otherwise we expect a description of the location. Wphase is
            # then run and results pushed pushed to the messaging system
            # and written to disk.

            # depth only needed for BoM XML
            try: depth = float(self.commandline().optionString("depth"))
            except Exception: depth = 0.

            try:
                self.eqinfo = {
                    'lat' : float(self.commandline().optionString("lat")),
                    'lon' : float(self.commandline().optionString("lon")),
                    'dep' : depth,
                    'time': UTCDateTime(self.commandline().optionString("time"))}
            except Exception:
                Logging.error('You must provide either lat/lon/time or a JSON payload')
                return False

            getter('sourcezone')
            getter('magtype', 'mag_type')
            getter('magvalue', 'mag_value', conv=float)
            getter('outputs', 'output')
            getter('server')
            getter('networks')
            getter('region')
            getter('evid')
            getter('resultid')
            getter('notificationemail')
            getter('emailmethod', 'email_method')
            getter('fromemail')
            getter('emailawsregion', 'email_aws_region', 'us-west-2')
            getter('writeS3', 'write_s3', False, lambda x: True)
            getter('bucketname', 'bucket_name')
            getter('agency')
            getter('make-maps', 'make_maps', False, lambda x: True)
            getter('overwrite', 'overwrite', False, lambda x: True)

            getter('smtp-server', 'smtp_server')
            getter('smtp-port', 'smtp_port')
            getter('smtp-user', 'smtp_user')
            getter('smtp-password', 'smtp_password')
            getter('smtp-ssl', 'smtp_ssl')
            getter('smtp-tls', 'smtp_tls')

            if self.evid is not None:
                self.output = os.path.join(self.output, self.evid)

            if self.resultid is not None:
                self.output = os.path.join(self.output, self.resultid)

            if self.notificationemail is not None:
                self.notificationemail = str(self.notificationemail).split(',')
                if not self.write_s3:
                    msg = 'requsted to send notification but not write to S3: will try and write to S3'
                    Logging.warning(msg)
                    self.write_s3 = True

                if not self.fromemail:
                    self.fromemail = self.notificationemail[0]

            if self.write_s3 and (
                    write_to_s3 is None or \
                    self.evid is None or \
                    self.bucket_name is None):
                Logging.error('attempt to write to s3, but no evid provided.')
                return False

            if self.notificationemail is not None and (
                    write_to_s3 is None or \
                    self.mag_type is None or \
                    self.mag_value is None or \
                    self.evid is None or \
                    self.resultid is None
                    ):
                Logging.error('cannot send email.')
                return False

        return True



    def init(self):
        """
        SC3 specific method.

        Returning False means that we do not enter the SeisComP3 run loop.
        """

        if Application.init(self) == False:
            return False

        item = None
        res = None
        parser = WPhaseParser(Logging.info)

        Logging.enableConsoleLogging(Logging.getGlobalChannel("error"))
        wphase_failed = False

        if self.evid is not None:
            self.eqinfo['id'] = self.evid
        if self.filename is None:
            try:
                wphase.logger.info("Starting W-Phase.")
                res = runwphase(
                    output_dir=self.output,
                    server=self.server,
                    eqinfo=self.eqinfo,
                    networks=self.networks,
                    make_maps=self.make_maps,
                    output_dir_can_exist=self.overwrite)
            except Exception:
                from traceback import format_exc
                Logging.error('failed to run wphase: {}'.format(format_exc()))
                wphase_failed = True
            else:
                if self.evid is not None:
                    try:
                        # TODO: Should this be done in runwphase?
                        res[settings.WPHASE_EVENT_KEY]['id'] = self.evid
                        with open(os.path.join(
                                self.output,
                                settings.WPHASE_OUTPUT_FILE_NAME), 'w') as of:
                            json.dump(res, of)
                    except Exception as e:
                        # not sure how we would get here, but we just don't want
                        # to stop the rest of processing
                        Logging.error('failed to add event id to event: {}'.format(e))

                try:
                    try: res_dict = res.as_dict()
                    except Exception: res_dict = res
                    item = parser.read(json_data=res_dict)
                except Exception as e:
                    Logging.error('failed to parse event JSON for SC3: {}'.format(e))
                    item = None

        else:
            try:
                item = parser.read(filename=self.filename)
            except Exception as e:
                Logging.error('failed parse event JSON for SC3: {}'.format(e))

        if item is not None:
            try:
                objs = createAndSendObjects(item, self.connection(),
                                            evid=self.evid, logging=Logging,
                                            agency=self.agency)
            except Exception as e:
                Logging.error('failed create objects for SC3: {}'.format(e))
            else:
                try:
                    # write output to file
                    filename = os.path.join(self.output, 'sc3.xml')
                    writeSCML(filename, objs)
                    Logging.info("Stored results in SCML as {}".format(filename))
                except Exception as e:
                    Logging.error('failed write SCML to file: {}'.format(e))


        if self.write_s3:
            # We always want to write to S3, even in the event of failure as
            # the JSON output may explain the failure. We also want to do it
            # at the very end since we keep the sc3 log file.

            try:
                write_to_s3(
                    self.output,
                    self.bucket_name,
                    self.evid,
                    self.resultid,
                    [(LOG_FILE_NAME, 'sc3.log')],
                    Logging.error)
            except Exception as e:
                Logging.error('failed write to S3: {}'.format(e))
            finally:
                # since we may do more runs, remove the log so we don't
                # contaminate later runs with info for this run
                try:
                    os.remove(LOG_FILE_NAME)
                except OSError:
                    # ... don't log this
                    pass

        if self.notificationemail:
            # must be done after writing to S3
            success = res is not None and 'MomentTensor' in res
            subject, body = self.createEmail(event_id=self.evid,
                                             result_id=self.resultid,
                                             result_dict=res,
                                             call_succeeded=success,
                                             )
            send_email(recipients=self.notificationemail,
                       subject=subject,
                       message=body,
                       from_email=self.fromemail,
                       method=self.email_method,

                       email_aws_region=self.email_aws_region,

                       server=self.smtp_server,
                       port=self.smtp_port,
                       user=self.smtp_user,
                       password=self.smtp_password,
                       ssl=self.smtp_ssl,
                       tls=self.smtp_tls,
                       )

        sys.exit(1 if wphase_failed else 0)

    def createEmail(self, event_id, result_id, result_dict, call_succeeded):
        """Create email subject and body to notify result.

        Override this in a subclass to produce a custom email."""
        subject = 'W-Phase succeeded' if call_succeeded else 'W-Phase failed'
        body = "W-Phase inversion {} for {} results: {}".format(result_id, event_id, result_dict)
        return subject, body


def main():
    argv = sys.argv
    app = WPhase(len(argv), argv)
    sys.exit(app())


if __name__ == "__main__":
    main()
